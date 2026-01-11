/**
 * Mol* Plugin Module
 *
 * This module provides full access to Mol* internal APIs via ES module imports.
 * Key APIs: StateTransforms, StructureElement.Bundle for creating sphere overlays.
 */

// Import Mol* CSS
import 'molstar/lib/mol-plugin-ui/skin/light.scss';

// Core imports
import { createPluginUI } from 'molstar/lib/mol-plugin-ui';
import { renderReact18 } from 'molstar/lib/mol-plugin-ui/react18';
import { DefaultPluginUISpec } from 'molstar/lib/mol-plugin-ui/spec';
import { PluginUIContext } from 'molstar/lib/mol-plugin-ui/context';
import { PluginContext } from 'molstar/lib/mol-plugin/context';
import { StateTransforms } from 'molstar/lib/mol-plugin-state/transforms';
import { StructureElement } from 'molstar/lib/mol-model/structure';
import { Color } from 'molstar/lib/mol-util/color';

// Export for use in main app
export { PluginContext, PluginUIContext, StateTransforms, StructureElement, Color };

// Global plugin instance
let plugin: PluginUIContext | null = null;

// Background color (cream/beige: rgb(228, 224, 221))
const BG_COLOR = 0xE4E0DD;

// Sphere halo tag for state tree queries
const SPHERE_HALO_TAG = 'residue-highlight-sphere';

// Track overlay tag prefix for multiple sphere overlays
const TRACK_OVERLAY_TAG_PREFIX = 'track-overlay-';

// Global ref for the sphere state node
let sphereHaloRef: string | null = null;

// Track overlay refs: Map of trackId -> array of state refs
const trackOverlayRefs: Map<string, string[]> = new Map();

// Debug flag
const DEBUG_HALO = true;

/**
 * Initialize the Mol* plugin with minimal UI
 */
export async function initMolstarPlugin(container: HTMLElement): Promise<PluginUIContext> {
    if (plugin) {
        await plugin.dispose();
    }

    // Create plugin with customized spec to hide UI elements
    const spec = DefaultPluginUISpec();

    // Hide the animation/controls buttons
    spec.layout = {
        initial: {
            isExpanded: false,
            showControls: false,
            regionState: {
                left: 'collapsed',
                right: 'hidden',
                top: 'hidden',
                bottom: 'hidden',
            }
        }
    };

    // Disable animation controls
    spec.components = {
        ...spec.components,
        remoteState: 'none'
    };

    plugin = await createPluginUI({
        target: container,
        spec: spec,
        render: renderReact18
    });

    // Wait for canvas to be ready
    await new Promise(resolve => setTimeout(resolve, 100));

    // Set background color via Canvas3D settings
    if (plugin.canvas3d) {
        plugin.canvas3d.setProps({
            renderer: {
                backgroundColor: Color(BG_COLOR),
            }
        });
    }

    // Expose for debugging
    (window as any).molstarPlugin = plugin;
    (window as any).Molstar = { StateTransforms, StructureElement };

    console.log('[Mol*] Plugin initialized with ES module imports.');
    console.log('[Mol*] StateTransforms:', StateTransforms);
    console.log('[Mol*] StructureElement.Bundle:', StructureElement.Bundle);

    return plugin;
}

/**
 * Get the current plugin instance
 */
export function getPlugin(): PluginUIContext | null {
    return plugin;
}

/**
 * Load a structure file into the Mol* viewer
 */
export async function loadStructure(url: string, format: string = 'mmcif'): Promise<void> {
    if (!plugin) throw new Error('Plugin not initialized');

    const data = await plugin.builders.data.download(
        { url, isBinary: false },
        { state: { isGhost: true } }
    );

    const trajectory = await plugin.builders.structure.parseTrajectory(data, format as any);

    await plugin.builders.structure.hierarchy.applyPreset(trajectory, 'default', {
        structure: {
            name: 'model',
            params: {}
        },
        showUnitcell: false,
        representationPreset: 'auto'
    });
}

/**
 * Remove existing sphere halo from the state tree
 */
export async function removeHalo(): Promise<void> {
    if (!plugin) {
        sphereHaloRef = null;
        return;
    }

    try {
        const state = plugin.state.data;

        // Method 1: Delete by stored ref (for state tree components)
        if (sphereHaloRef) {
            const cell = state.cells.get(sphereHaloRef);
            if (cell) {
                await plugin.build().delete(sphereHaloRef).commit();
                if (DEBUG_HALO) {
                    console.log('[HALO] Sphere halo removed by ref:', sphereHaloRef);
                }
            }
        }

        // Method 2: Find by tag and delete (most reliable for state tree nodes)
        const toDelete: string[] = [];
        state.cells.forEach((cell, ref) => {
            if (cell.transform?.tags?.includes(SPHERE_HALO_TAG)) {
                toDelete.push(ref);
            }
        });

        if (toDelete.length > 0) {
            const update = plugin.build();
            for (const ref of toDelete) {
                update.delete(ref);
            }
            await update.commit();
            if (DEBUG_HALO) {
                console.log('[HALO] Removed', toDelete.length, 'sphere halo(s) by tag');
            }
        }

        // Method 3: Clear structure selection
        if (plugin.managers?.structure?.selection) {
            plugin.managers.structure.selection.clear();
        }

        sphereHaloRef = null;

    } catch (e) {
        if (DEBUG_HALO) {
            console.log('[HALO] Error in removeHalo:', (e as Error).message);
        }
        sphereHaloRef = null;
    }
}

/**
 * Extract a single-atom loci (CA preferred, or first atom fallback) from a residue loci.
 */
export function extractCentralAtomLoci(residueLoci: any): any {
    if (!residueLoci || !residueLoci.elements?.length) {
        return residueLoci;
    }

    try {
        for (const e of residueLoci.elements) {
            const unit = e.unit;
            if (!unit?.model) continue;

            const { atomicHierarchy } = unit.model;
            const { label_atom_id } = atomicHierarchy.atoms;
            const originalIndices = e.indices;

            // Extract indices array
            let indicesArray: number[] = [];

            if (originalIndices && typeof originalIndices.length === 'number') {
                for (let i = 0; i < originalIndices.length; i++) {
                    indicesArray.push(originalIndices[i]);
                }
            } else if (typeof originalIndices === 'number') {
                // Encoded interval format (legacy)
                const min = originalIndices >>> 16;
                const max = originalIndices & 0xFFFF;
                for (let i = min; i < max; i++) {
                    indicesArray.push(i);
                }
            } else {
                continue;
            }

            if (indicesArray.length === 0) continue;

            // Search for CA atom first
            let caIdx: number | null = null;
            let c4PrimeIdx: number | null = null;
            let firstIdx = indicesArray[0];

            for (const idx of indicesArray) {
                const atomIdx = unit.elements[idx];
                if (atomIdx === undefined) continue;

                const atomName = label_atom_id.value(atomIdx);

                if (atomName === 'CA') {
                    caIdx = idx;
                    break;
                }
                if (atomName === "C4'" && c4PrimeIdx === null) {
                    c4PrimeIdx = idx;
                }
            }

            const targetIdx = caIdx !== null ? caIdx : (c4PrimeIdx !== null ? c4PrimeIdx : firstIdx);

            // Create a proper SortedArray with single element
            let singleAtomIndices: Int32Array;

            const IndicesConstructor = originalIndices?.constructor as any;
            if (IndicesConstructor?.ofSingleton) {
                singleAtomIndices = IndicesConstructor.ofSingleton(targetIdx);
            } else {
                singleAtomIndices = new Int32Array([targetIdx]);
            }

            // Construct the new single-atom loci
            const singleAtomLoci = {
                kind: 'element-loci',
                structure: residueLoci.structure,
                elements: [{
                    unit: unit,
                    indices: singleAtomIndices
                }]
            };

            if (DEBUG_HALO) {
                const atomIdx = unit.elements[targetIdx];
                const atomName = label_atom_id.value(atomIdx);
                console.log(`[HALO] extractCentralAtomLoci: selected atom ${atomName} at index ${targetIdx}`);
            }

            return singleAtomLoci;
        }

        return residueLoci;

    } catch (err) {
        console.log('[HALO] extractCentralAtomLoci error:', (err as Error).message);
        return residueLoci;
    }
}

/**
 * Create a sphere halo around the clicked residue using single-atom spacefill
 * This approach creates a large spacefill representation on just the CA atom
 * to create a "halo sphere" effect.
 */
export async function createSphereHalo(
    loci: any,
    residueNum: number,
    sizeFactor: number = 6.0,
    alpha: number = 0.4
): Promise<boolean> {
    console.log('[HALO] createSphereHalo called, loci:', loci, 'residueNum:', residueNum);

    if (!loci || loci.elements?.length === 0) {
        console.warn('[HALO] empty loci, abort halo');
        return false;
    }

    if (!plugin) {
        console.warn('[HALO] plugin not available');
        return false;
    }

    try {
        // 1. Remove existing halo first
        console.log('[HALO] removing existing halo');
        await removeHalo();

        // 2. Extract single-atom loci (CA preferred, fallback to first atom)
        console.log('[HALO] extracting central atom loci');
        const singleAtomLoci = extractCentralAtomLoci(loci);
        if (!singleAtomLoci || !singleAtomLoci.elements?.length) {
            console.warn('[HALO] could not extract single-atom loci');
            return false;
        }
        console.log('[HALO] singleAtomLoci:', singleAtomLoci);

        // 3. Get the structure from the loci
        const structure = singleAtomLoci.structure;
        if (!structure) {
            console.warn('[HALO] no structure in loci');
            return false;
        }
        console.log('[HALO] structure found');

        // 4. Find the structure cell in the state tree
        let structureRef: string | null = null;
        plugin.state.data.cells.forEach((cell, ref) => {
            if (cell.obj?.data === structure) {
                structureRef = ref;
            }
        });

        if (!structureRef) {
            console.warn('[HALO] could not find structure in state tree');
            return false;
        }
        console.log('[HALO] structureRef:', structureRef);

        // 5. Create spacefill halo using StructureElement.Bundle
        const Bundle = StructureElement.Bundle;

        console.log('[HALO] StateTransforms:', StateTransforms);
        console.log('[HALO] StructureElement:', StructureElement);
        console.log('[HALO] Bundle:', Bundle);

        // Yellow halo color
        const HALO_COLOR = Color(0xFFE632);

        if (Bundle?.fromLoci && StateTransforms?.Model?.StructureSelectionFromBundle) {
            console.log('[HALO] building selection from loci via Bundle');
            const bundle = Bundle.fromLoci(singleAtomLoci);
            console.log('[HALO] bundle created:', bundle);

            // Build representation params with transparency
            const reprParams = {
                type: {
                    name: 'spacefill' as const,
                    params: {
                        sizeFactor: sizeFactor,
                        alpha: alpha
                    }
                },
                colorTheme: { name: 'uniform' as const, params: { value: HALO_COLOR } },
                sizeTheme: { name: 'uniform' as const, params: { value: sizeFactor } }
            };
            console.log('[HALO] reprParams:', reprParams);

            const result = await plugin.build()
                .to(structureRef)
                .apply(StateTransforms.Model.StructureSelectionFromBundle, {
                    bundle: bundle,
                    label: `Residue ${residueNum} Halo`
                }, { tags: [SPHERE_HALO_TAG] })
                .apply(StateTransforms.Representation.StructureRepresentation3D, reprParams, { tags: [SPHERE_HALO_TAG] })
                .commit();

            sphereHaloRef = result?.ref || 'spacefill-halo-' + residueNum;
            console.log('[HALO] SUCCESS via StructureSelectionFromBundle, ref:', sphereHaloRef);
            return true;
        }

        console.warn('[HALO] Bundle API not available - no suitable API found');
        return false;

    } catch (e) {
        console.error('[HALO] Failed to create sphere halo:', (e as Error).message, (e as Error).stack);
        return false;
    }
}

/**
 * Reset camera to initial view
 */
export function resetCamera(): void {
    if (!plugin?.canvas3d?.camera) return;
    plugin.canvas3d.requestCameraReset();
}

/**
 * Subscribe to click events on residues
 */
export function subscribeToClicks(callback: (loci: any, residueNum: number) => void): void {
    if (!plugin) return;

    plugin.behaviors.interaction.click.subscribe(({ current, button }) => {
        if (button !== 0) return; // Left click only

        const loci = current.loci;
        if (!loci || loci.kind !== 'element-loci' || !loci.elements?.length) return;

        // Extract residue number from loci
        try {
            const element = loci.elements[0];
            const unit = element.unit;
            if (!unit?.model) return;

            const { residues, residueAtomSegments } = unit.model.atomicHierarchy;
            const { label_seq_id } = residues;

            // Get first atom index from the loci
            let atomIdx: number;
            const indices = element.indices as any;
            if (typeof indices === 'number') {
                atomIdx = indices >>> 16;
            } else if (indices && (indices.length > 0 || indices[0] !== undefined)) {
                atomIdx = indices[0];
            } else {
                return;
            }

            const unitAtomIdx = unit.elements[atomIdx];
            const residueIndex = residueAtomSegments.index[unitAtomIdx];
            const residueNum = label_seq_id.value(residueIndex);

            callback(loci, residueNum);
        } catch (e) {
            console.error('[Mol*] Error extracting residue from click:', e);
        }
    });
}

/**
 * Disable default Mol* behaviors that interfere with custom handling
 */
export function disableDefaultBehaviors(): void {
    if (!plugin) return;

    try {
        // Disable marking/highlighting via Canvas3D settings
        if (plugin.canvas3d?.setProps) {
            plugin.canvas3d.setProps({
                marking: {
                    enabled: false,
                },
                renderer: {
                    ...plugin.canvas3d.props?.renderer,
                    selectMark: false,
                    highlightStrength: 0,
                    selectStrength: 0,
                },
                highlightColor: Color.fromRgb(0, 0, 0),
                selectColor: Color.fromRgb(0, 0, 0),
            } as any);
        }
    } catch (e) {
        console.warn('[Mol*] Could not disable default behaviors:', e);
    }
}

/**
 * Create track overlay spheres for multiple residues
 * Used for visualizing variant tracks (ClinVar, Training, Constraint) on the 3D structure
 *
 * @param residueNumbers Array of residue numbers to highlight
 * @param color Hex color value for the spheres
 * @param trackId Unique identifier for this track overlay
 * @param sizeFactor Size multiplier for spheres (default 3.0)
 * @param alpha Transparency (default 0.7)
 * @returns Array of state refs for the created spheres
 */
export async function createTrackOverlaySpheres(
    residueNumbers: number[],
    color: number,
    trackId: string,
    sizeFactor: number = 3.0,
    alpha: number = 0.7
): Promise<string[]> {
    if (!plugin) {
        console.warn('[OVERLAY] Plugin not available');
        return [];
    }

    if (!residueNumbers || residueNumbers.length === 0) {
        console.warn('[OVERLAY] No residue numbers provided');
        return [];
    }

    const tag = TRACK_OVERLAY_TAG_PREFIX + trackId;
    const sphereColor = Color(color);
    const refs: string[] = [];

    try {
        // Remove any existing overlay for this track
        await removeTrackOverlaySpheres(trackId);

        // Find the structure in the state tree
        let structureRef: string | null = null;
        let structure: any = null;

        plugin.state.data.cells.forEach((cell, ref) => {
            if (cell.obj?.type?.name === 'Structure' && cell.obj?.data) {
                structureRef = ref;
                structure = cell.obj.data;
            }
        });

        if (!structureRef || !structure) {
            console.warn('[OVERLAY] No structure found in state tree');
            return [];
        }

        console.log(`[OVERLAY] Creating spheres for ${residueNumbers.length} residues, track: ${trackId}`);

        // Build selection expression for all residues
        // Use residue-level selection with spacefill representation
        const Bundle = StructureElement.Bundle;

        if (!Bundle?.fromLoci || !StateTransforms?.Model?.StructureSelectionFromBundle) {
            console.warn('[OVERLAY] Required APIs not available');
            return [];
        }

        // Create a single combined loci for all residues
        // This is more efficient than creating individual spheres
        const model = structure.models?.[0];
        if (!model) {
            console.warn('[OVERLAY] No model found in structure');
            return [];
        }

        // Build the update with representation for selected residues
        const reprParams = {
            type: {
                name: 'spacefill' as const,
                params: {
                    sizeFactor: sizeFactor,
                    alpha: alpha
                }
            },
            colorTheme: { name: 'uniform' as const, params: { value: sphereColor } },
            sizeTheme: { name: 'uniform' as const, params: { value: sizeFactor } }
        };

        // Create selection query for the residues
        // Use MolScript-style residue selection
        const { MolScriptBuilder: MS } = await import('molstar/lib/mol-script/language/builder');
        const { compile } = await import('molstar/lib/mol-script/script');
        const { Script } = await import('molstar/lib/mol-script/script');

        // Build residue selection expression
        const residueSet = new Set(residueNumbers);
        const residueExpr = MS.struct.generator.atomGroups({
            'residue-test': MS.core.set.has([
                MS.set(...Array.from(residueSet)),
                MS.ammp('auth_seq_id')
            ])
        });

        const result = await plugin.build()
            .to(structureRef)
            .apply(StateTransforms.Model.StructureSelectionFromExpression, {
                expression: residueExpr,
                label: `Track Overlay: ${trackId}`
            }, { tags: [tag] })
            .apply(StateTransforms.Representation.StructureRepresentation3D, reprParams, { tags: [tag] })
            .commit();

        if (result?.ref) {
            refs.push(result.ref);
        }

        // Store refs for this track
        trackOverlayRefs.set(trackId, refs);

        console.log(`[OVERLAY] Created overlay for track ${trackId}, refs:`, refs);
        return refs;

    } catch (e) {
        console.error('[OVERLAY] Failed to create track overlay:', (e as Error).message, (e as Error).stack);
        return [];
    }
}

/**
 * Remove track overlay spheres for a specific track
 *
 * @param trackId The track identifier to remove overlays for
 */
export async function removeTrackOverlaySpheres(trackId: string): Promise<void> {
    if (!plugin) {
        trackOverlayRefs.delete(trackId);
        return;
    }

    const tag = TRACK_OVERLAY_TAG_PREFIX + trackId;

    try {
        const state = plugin.state.data;
        const toDelete: string[] = [];

        // Find all cells with this track's tag
        state.cells.forEach((cell, ref) => {
            if (cell.transform?.tags?.includes(tag)) {
                toDelete.push(ref);
            }
        });

        if (toDelete.length > 0) {
            const update = plugin.build();
            for (const ref of toDelete) {
                update.delete(ref);
            }
            await update.commit();
            console.log(`[OVERLAY] Removed ${toDelete.length} overlay(s) for track ${trackId}`);
        }

        trackOverlayRefs.delete(trackId);

    } catch (e) {
        console.warn('[OVERLAY] Error removing track overlay:', (e as Error).message);
        trackOverlayRefs.delete(trackId);
    }
}

/**
 * Clear all track overlays from the structure
 */
export async function clearAllTrackOverlays(): Promise<void> {
    if (!plugin) {
        trackOverlayRefs.clear();
        return;
    }

    try {
        const state = plugin.state.data;
        const toDelete: string[] = [];

        // Find all cells with track overlay tags
        state.cells.forEach((cell, ref) => {
            const tags = cell.transform?.tags || [];
            for (const t of tags) {
                if (typeof t === 'string' && t.startsWith(TRACK_OVERLAY_TAG_PREFIX)) {
                    toDelete.push(ref);
                    break;
                }
            }
        });

        if (toDelete.length > 0) {
            const update = plugin.build();
            for (const ref of toDelete) {
                update.delete(ref);
            }
            await update.commit();
            console.log(`[OVERLAY] Cleared ${toDelete.length} track overlay(s)`);
        }

        trackOverlayRefs.clear();

    } catch (e) {
        console.warn('[OVERLAY] Error clearing track overlays:', (e as Error).message);
        trackOverlayRefs.clear();
    }
}

/**
 * Check if a track currently has an overlay
 *
 * @param trackId The track identifier to check
 * @returns true if the track has an active overlay
 */
export function hasTrackOverlay(trackId: string): boolean {
    return trackOverlayRefs.has(trackId) && (trackOverlayRefs.get(trackId)?.length ?? 0) > 0;
}

/**
 * Get all currently overlayed track IDs
 *
 * @returns Array of track IDs with active overlays
 */
export function getOverlayedTracks(): string[] {
    return Array.from(trackOverlayRefs.keys());
}
