# Variant Track Protein Overlay Feature

## Feature Overview

This feature enables users to overlay "variant" tracks from the D3 genomic track viewer onto the 3D protein structure in Mol*, displaying positions with data as colored spheres on the protein surface.

### Requirements Summary

| Requirement | Specification |
|-------------|---------------|
| Max simultaneous overlays | 2 tracks |
| Track 1 color | `#0066CC` (blue) |
| Track 2 color | `#00CED1` (dark cyan) |
| Eligible track types | ClinVar, Training Labels, Constraint/Allele |
| Position filter | Non-zero values only |
| Filter reactivity | Spheres update when track filters change |
| Yellow highlight | Must coexist with track spheres |
| UI control location | Button below flip button on track colorbar |

---

## Eligible Tracks

The overlay feature is only available for "variant" type tracks (discrete/categorical data), not continuous tracks.

### Detection Logic

```javascript
function isEligibleForSphereOverlay(trackId) {
    return isClinVarStackedTrack(trackId) ||      // clinvar.clinvar_label_list
           isTrainingTrack(trackId) ||             // training.train_counts.*
           isConstraintStackedTrack(trackId);      // Constraint, Core, Complete
}
```

### Track Data Structures

| Track Type | Field ID | Data Structure | Non-Zero Check |
|------------|----------|----------------|----------------|
| ClinVar Labels | `clinvar.clinvar_label_list` | `Array<string>` | `labels.filter(clinvarLabelPassesFilter).length > 0` |
| Training Labelled | `training.train_counts.labelled` | `number` | `count > 0` |
| Training Unlabelled | `training.train_counts.unlabelled` | `number` | `count > 0` |
| Training Labelled HQ | `training.train_counts.labelled_high_qual` | `number` | `count > 0` |
| Training Unlabelled HQ | `training.train_counts.unlabelled_high_qual` | `number` | `count > 0` |
| Constraint | `Constraint` | `Array<{_0: allele, _1: pred}>` | `variants.length > 0` |
| Core | `Core` | `Array<{_0: allele, _1: pred}>` | `variants.length > 0` |
| Complete | `Complete` | `Array<{_0: allele, _1: pred}>` | `variants.length > 0` |

---

## Architecture

### Data Flow

```
┌─────────────────────────────────────────────────────────────────────────┐
│                     USER CLICKS OVERLAY BUTTON                          │
│                        (on eligible track)                              │
└──────────────────────────────┬──────────────────────────────────────────┘
                               │
                               ▼
┌─────────────────────────────────────────────────────────────────────────┐
│                     CHECK OVERLAY CAPACITY                              │
│                                                                         │
│  canAddTrackOverlay() → returns true if < 2 tracks overlayed           │
│  If false, button is disabled                                          │
└──────────────────────────────┬──────────────────────────────────────────┘
                               │
                               ▼
┌─────────────────────────────────────────────────────────────────────────┐
│                     ASSIGN COLOR SLOT                                   │
│                                                                         │
│  Slot 1 (empty) → #0066CC (blue)                                       │
│  Slot 2 (if slot 1 taken) → #00CED1 (cyan)                             │
└──────────────────────────────┬──────────────────────────────────────────┘
                               │
                               ▼
┌─────────────────────────────────────────────────────────────────────────┐
│                   EXTRACT RESIDUES WITH DATA                            │
│                                                                         │
│  getResiduesWithTrackData(trackId):                                    │
│    1. Iterate positionsData (current visible window)                   │
│    2. For each position, check if has non-zero value for trackId       │
│    3. Apply track-specific filters (e.g., clinvarLabelFilter)          │
│    4. Map genomic position → protein residue (aa_pos column)           │
│    5. Deduplicate (multiple codons → same amino acid)                  │
│    6. Return Array<residueNumber>                                      │
└──────────────────────────────┬──────────────────────────────────────────┘
                               │
                               ▼
┌─────────────────────────────────────────────────────────────────────────┐
│                   CREATE SPHERES IN MOL*                                │
│                                                                         │
│  createTrackOverlaySpheres(residues, color, trackId):                  │
│    1. Find structure ref in Mol* state tree                            │
│    2. For each residue, create selection expression                    │
│    3. Apply spacefill representation with:                             │
│       - sizeFactor: 3.0 (smaller than yellow halo)                     │
│       - alpha: 0.7                                                      │
│       - uniform color from slot                                        │
│    4. Tag all with: track-overlay-{trackId}                            │
│    5. Batch commit for performance                                     │
└─────────────────────────────────────────────────────────────────────────┘
```

### Filter Change Flow

```
┌─────────────────────────────────────────────────────────────────────────┐
│                  USER CHANGES CLINVAR FILTER                            │
│                (toggles Pathogenic/Benign/VUS/etc.)                     │
└──────────────────────────────┬──────────────────────────────────────────┘
                               │
                               ▼
┌─────────────────────────────────────────────────────────────────────────┐
│              CHECK IF CLINVAR TRACK IS OVERLAYED                        │
│                                                                         │
│  if (isTrackOverlayed('clinvar.clinvar_label_list'))                   │
└──────────────────────────────┬──────────────────────────────────────────┘
                               │
                               ▼
┌─────────────────────────────────────────────────────────────────────────┐
│                   UPDATE TRACK SPHERES                                  │
│                                                                         │
│  updateTrackOverlaySpheres('clinvar.clinvar_label_list'):              │
│    1. removeTrackOverlaySpheres(trackId)                               │
│    2. getResiduesWithTrackData(trackId) // uses new filter             │
│    3. createTrackOverlaySpheres(residues, color, trackId)              │
└─────────────────────────────────────────────────────────────────────────┘
```

---

## Implementation Details

### 1. State Management

**Location:** `gosling_mvp/frontend/d3_viewer.html` (near line 1800, after existing state variables)

```javascript
// ============================================
// TRACK OVERLAY STATE
// ============================================

/**
 * State for track sphere overlays on 3D structure.
 * Maximum 2 tracks can be overlayed simultaneously.
 */
let trackOverlayState = {
    // Array of currently overlayed tracks: [{trackId, colorSlot}]
    // colorSlot: 0 = #0066CC, 1 = #00CED1
    overlayedTracks: [],

    // Cache of created sphere refs for each track (for cleanup)
    sphereRefs: new Map()  // Map<trackId, string[]>
};

// Fixed overlay colors (slot index → hex color)
const TRACK_OVERLAY_COLORS = [
    0x0066CC,  // Slot 0: Blue
    0x00CED1   // Slot 1: Dark Cyan
];

// CSS versions for UI
const TRACK_OVERLAY_COLORS_CSS = [
    '#0066CC',
    '#00CED1'
];

/**
 * Check if a track type is eligible for sphere overlay.
 * Only "variant" tracks (discrete data) are eligible.
 */
function isEligibleForSphereOverlay(trackId) {
    return isClinVarStackedTrack(trackId) ||
           isTrainingTrack(trackId) ||
           isConstraintStackedTrack(trackId);
}

/**
 * Check if a track is currently overlayed on the 3D structure.
 */
function isTrackOverlayed(trackId) {
    return trackOverlayState.overlayedTracks.some(t => t.trackId === trackId);
}

/**
 * Check if we can add another track overlay (max 2).
 */
function canAddTrackOverlay() {
    return trackOverlayState.overlayedTracks.length < 2;
}

/**
 * Get the next available color slot (0 or 1).
 * Returns null if no slots available.
 */
function getNextColorSlot() {
    if (trackOverlayState.overlayedTracks.length === 0) return 0;
    if (trackOverlayState.overlayedTracks.length === 1) {
        // Return the slot NOT used by the first track
        return trackOverlayState.overlayedTracks[0].colorSlot === 0 ? 1 : 0;
    }
    return null;  // Both slots taken
}

/**
 * Get the color for a track's overlay (hex number for Mol*).
 * Returns null if track is not overlayed.
 */
function getTrackOverlayColor(trackId) {
    const track = trackOverlayState.overlayedTracks.find(t => t.trackId === trackId);
    if (!track) return null;
    return TRACK_OVERLAY_COLORS[track.colorSlot];
}

/**
 * Get the CSS color for a track's overlay.
 */
function getTrackOverlayColorCSS(trackId) {
    const track = trackOverlayState.overlayedTracks.find(t => t.trackId === trackId);
    if (!track) return null;
    return TRACK_OVERLAY_COLORS_CSS[track.colorSlot];
}
```

### 2. Mol* Plugin Functions

**Location:** `gosling_mvp/frontend/src/molstar-plugin.ts`

```typescript
// ============================================
// TRACK OVERLAY SPHERE MANAGEMENT
// ============================================

// Tag prefix for track overlay spheres (distinct from click highlight)
const TRACK_OVERLAY_TAG_PREFIX = 'track-overlay-';

/**
 * Create multiple spheres on the protein structure for residues with track data.
 * Uses a single batch commit for performance.
 *
 * @param residueNumbers - Array of protein residue numbers to highlight
 * @param color - Hex color (e.g., 0x0066CC)
 * @param trackId - Track identifier for tagging (enables cleanup)
 * @param sizeFactor - Sphere size (default 3.0, smaller than yellow halo)
 * @param alpha - Transparency (default 0.7)
 * @returns Array of created state refs
 */
export async function createTrackOverlaySpheres(
    residueNumbers: number[],
    color: number,
    trackId: string,
    sizeFactor: number = 3.0,
    alpha: number = 0.7
): Promise<string[]> {
    if (!plugin) {
        console.warn('[TRACK-OVERLAY] Plugin not initialized');
        return [];
    }

    if (residueNumbers.length === 0) {
        console.log('[TRACK-OVERLAY] No residues to highlight');
        return [];
    }

    const tag = `${TRACK_OVERLAY_TAG_PREFIX}${trackId}`;
    console.log(`[TRACK-OVERLAY] Creating ${residueNumbers.length} spheres for ${trackId}`);

    try {
        // Find the structure in state tree
        let structureRef: string | null = null;
        let structureData: any = null;

        plugin.state.data.cells.forEach((cell, ref) => {
            if (cell.obj?.type?.name === 'Structure' && !structureRef) {
                structureRef = ref;
                structureData = cell.obj.data;
            }
        });

        if (!structureRef || !structureData) {
            console.warn('[TRACK-OVERLAY] No structure found in state tree');
            return [];
        }

        // Build selection expression for all residues
        // Using MS.struct.generator.atomGroups with auth_seq_id filter
        const MS = (window as any).molstar?.StructureSelectionQuery?.MS;

        // Create a selection for all specified residues at once
        // This is more efficient than creating individual spheres
        const residueSet = new Set(residueNumbers);

        const reprParams = {
            type: {
                name: 'spacefill' as const,
                params: {
                    sizeFactor: sizeFactor,
                    alpha: alpha
                }
            },
            colorTheme: {
                name: 'uniform' as const,
                params: { value: Color(color) }
            },
            sizeTheme: {
                name: 'uniform' as const,
                params: { value: sizeFactor }
            }
        };

        // Create selection and representation
        // Use StructureSelectionFromExpression for residue-based selection
        const result = await plugin.build()
            .to(structureRef)
            .apply(StateTransforms.Model.StructureSelectionFromExpression, {
                expression: MS.struct.generator.atomGroups({
                    'residue-test': MS.core.set.has([
                        MS.set(...residueNumbers),
                        MS.ammp('label_seq_id')
                    ]),
                    'atom-test': MS.core.rel.eq([
                        MS.ammp('label_atom_id'),
                        'CA'  // Only CA atoms for sphere centers
                    ])
                }),
                label: `Track Overlay: ${trackId}`
            }, { tags: [tag] })
            .apply(StateTransforms.Representation.StructureRepresentation3D, reprParams, { tags: [tag] })
            .commit();

        const refs = result?.ref ? [result.ref] : [];
        console.log(`[TRACK-OVERLAY] Created spheres with refs:`, refs);
        return refs;

    } catch (e) {
        console.error('[TRACK-OVERLAY] Failed to create spheres:', (e as Error).message);
        return [];
    }
}

/**
 * Remove all track overlay spheres for a specific track.
 *
 * @param trackId - Track identifier to remove spheres for
 */
export async function removeTrackOverlaySpheres(trackId: string): Promise<void> {
    if (!plugin) return;

    const tag = `${TRACK_OVERLAY_TAG_PREFIX}${trackId}`;
    console.log(`[TRACK-OVERLAY] Removing spheres for ${trackId}`);

    try {
        const state = plugin.state.data;
        const toDelete: string[] = [];

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
            console.log(`[TRACK-OVERLAY] Removed ${toDelete.length} nodes for ${trackId}`);
        }

    } catch (e) {
        console.error('[TRACK-OVERLAY] Failed to remove spheres:', (e as Error).message);
    }
}

/**
 * Remove all track overlay spheres from all tracks.
 */
export async function clearAllTrackOverlays(): Promise<void> {
    if (!plugin) return;

    console.log('[TRACK-OVERLAY] Clearing all track overlays');

    try {
        const state = plugin.state.data;
        const toDelete: string[] = [];

        state.cells.forEach((cell, ref) => {
            const tags = cell.transform?.tags || [];
            if (tags.some(t => t.startsWith(TRACK_OVERLAY_TAG_PREFIX))) {
                toDelete.push(ref);
            }
        });

        if (toDelete.length > 0) {
            const update = plugin.build();
            for (const ref of toDelete) {
                update.delete(ref);
            }
            await update.commit();
            console.log(`[TRACK-OVERLAY] Cleared ${toDelete.length} overlay nodes`);
        }

    } catch (e) {
        console.error('[TRACK-OVERLAY] Failed to clear overlays:', (e as Error).message);
    }
}
```

### 3. Data Extraction Function

**Location:** `gosling_mvp/frontend/d3_viewer.html` (after state management functions)

```javascript
/**
 * Get protein residue numbers that have non-zero data for a track.
 * Respects current filters (e.g., ClinVar significance filter).
 *
 * @param {string} trackId - The track field ID
 * @returns {number[]} - Array of unique protein residue numbers
 */
function getResiduesWithTrackData(trackId) {
    if (!currentStructureGene || !positionsData || positionsData.length === 0) {
        console.log('[OVERLAY] No structure or position data');
        return [];
    }

    const residuesWithData = new Set();

    for (const pos of positionsData) {
        // Skip positions without amino acid mapping
        // The aa_pos field should be populated by the backend
        const aaPos = pos.aa_pos || pos.protein_residue;
        if (!aaPos) continue;

        // Skip positions not in the current gene
        if (pos.gene_symbol?.toUpperCase() !== currentStructureGene.toUpperCase()) continue;

        let hasData = false;

        if (isClinVarStackedTrack(trackId)) {
            // ClinVar: check if any labels pass the current filter
            const labels = pos['clinvar.clinvar_label_list'];
            if (labels && Array.isArray(labels)) {
                const filteredLabels = labels.filter(label => clinvarLabelPassesFilter(label));
                hasData = filteredLabels.length > 0;
            }
        } else if (isTrainingTrack(trackId)) {
            // Training tracks: check if count > 0
            const count = pos[trackId];
            hasData = count !== null && count !== undefined && count > 0;
        } else if (isConstraintStackedTrack(trackId)) {
            // Constraint/Core/Complete: check if has variant predictions
            const variants = pos[trackId];
            hasData = variants && Array.isArray(variants) && variants.length > 0;
        }

        if (hasData) {
            residuesWithData.add(aaPos);
        }
    }

    const result = Array.from(residuesWithData).sort((a, b) => a - b);
    console.log(`[OVERLAY] Found ${result.length} residues with data for ${trackId}`);
    return result;
}
```

### 4. Overlay Control Functions

**Location:** `gosling_mvp/frontend/d3_viewer.html` (after getResiduesWithTrackData)

```javascript
/**
 * Add a track overlay to the 3D structure.
 * Creates spheres at residues with non-zero track values.
 *
 * @param {string} trackId - Track to overlay
 */
async function addTrackOverlay(trackId) {
    if (!plugin || !currentStructureGene) {
        console.warn('[OVERLAY] Cannot add overlay: no plugin or structure');
        return;
    }

    if (!canAddTrackOverlay()) {
        console.warn('[OVERLAY] Maximum 2 overlays reached');
        return;
    }

    if (isTrackOverlayed(trackId)) {
        console.warn('[OVERLAY] Track already overlayed:', trackId);
        return;
    }

    // Get next available color slot
    const colorSlot = getNextColorSlot();
    if (colorSlot === null) {
        console.warn('[OVERLAY] No color slot available');
        return;
    }

    // Add to state
    trackOverlayState.overlayedTracks.push({ trackId, colorSlot });

    // Get residues with data
    const residues = getResiduesWithTrackData(trackId);

    if (residues.length === 0) {
        console.log('[OVERLAY] No residues with data for track:', trackId);
        // Still keep in state so button shows as active
    } else {
        // Create spheres
        const color = TRACK_OVERLAY_COLORS[colorSlot];
        const refs = await createTrackOverlaySpheres(residues, color, trackId);
        trackOverlayState.sphereRefs.set(trackId, refs);
    }

    // Re-render to update button state
    renderChart();
}

/**
 * Remove a track overlay from the 3D structure.
 *
 * @param {string} trackId - Track to remove
 */
async function removeTrackOverlay(trackId) {
    if (!isTrackOverlayed(trackId)) {
        return;
    }

    // Remove spheres from Mol*
    await removeTrackOverlaySpheres(trackId);

    // Update state
    trackOverlayState.overlayedTracks = trackOverlayState.overlayedTracks
        .filter(t => t.trackId !== trackId);
    trackOverlayState.sphereRefs.delete(trackId);

    // Re-render to update button state
    renderChart();
}

/**
 * Update spheres for a track (called when filters change).
 *
 * @param {string} trackId - Track to update
 */
async function updateTrackOverlaySpheres(trackId) {
    if (!isTrackOverlayed(trackId)) return;

    // Get current color for this track
    const color = getTrackOverlayColor(trackId);
    if (color === null) return;

    // Remove existing spheres
    await removeTrackOverlaySpheres(trackId);

    // Get updated residue list (with new filter applied)
    const residues = getResiduesWithTrackData(trackId);

    // Create new spheres
    if (residues.length > 0) {
        const refs = await createTrackOverlaySpheres(residues, color, trackId);
        trackOverlayState.sphereRefs.set(trackId, refs);
    } else {
        trackOverlayState.sphereRefs.delete(trackId);
    }
}

/**
 * Clear all track overlays (e.g., when structure changes).
 */
async function clearTrackOverlays() {
    await clearAllTrackOverlays();
    trackOverlayState.overlayedTracks = [];
    trackOverlayState.sphereRefs.clear();
}
```

### 5. UI Button Implementation

**Location:** `gosling_mvp/frontend/d3_viewer.html` (in renderChart function, after flip button ~line 3061)

```javascript
// ============================================
// TRACK OVERLAY BUTTON
// ============================================

// Check if this track is eligible for sphere overlay
const isOverlayEligible = isEligibleForSphereOverlay(trackId);
const hasStructure = currentStructureGene && plugin;

if (isOverlayEligible && hasStructure) {
    const isOverlayed = isTrackOverlayed(trackId);
    const canAdd = canAddTrackOverlay();
    const overlayColorCSS = isOverlayed ? getTrackOverlayColorCSS(trackId) :
                            (canAdd ? TRACK_OVERLAY_COLORS_CSS[getNextColorSlot() || 0] : '#666');

    // Calculate button positions to center both buttons
    // If both buttons present: flip at -12, overlay at +12
    // If only flip button: flip at 0
    const hasFlipButton = true;  // Flip button always present for tracks with colorbars
    const buttonGroupOffset = hasFlipButton ? 0 : 0;
    const flipBtnY = -12;
    const overlayBtnY = 12;

    // Adjust flip button position if we're showing overlay button
    // (This should be done where flip button is created)

    // Create overlay button
    const overlayBtn = colorbarG.append('g')
        .attr('class', 'overlay-btn')
        .attr('transform', `translate(${colorbarWidth + 35}, ${colorbarHeight / 2 + overlayBtnY})`)
        .style('cursor', (isOverlayed || canAdd) ? 'pointer' : 'not-allowed')
        .style('opacity', (isOverlayed || canAdd) ? 1 : 0.4);

    // Button circle
    overlayBtn.append('circle')
        .attr('r', 10)
        .attr('fill', isOverlayed ? overlayColorCSS : '#333')
        .attr('stroke', overlayColorCSS)
        .attr('stroke-width', isOverlayed ? 2 : 1.5);

    // "3D" label
    overlayBtn.append('text')
        .attr('text-anchor', 'middle')
        .attr('dominant-baseline', 'middle')
        .attr('font-size', '8px')
        .attr('font-weight', 'bold')
        .attr('fill', isOverlayed ? '#fff' : overlayColorCSS)
        .attr('pointer-events', 'none')
        .text('3D');

    // Tooltip
    overlayBtn.append('title')
        .text(isOverlayed ?
              `Remove ${trackId} from 3D structure` :
              (canAdd ?
               `Overlay ${trackId} on 3D structure` :
               'Maximum 2 overlays allowed'));

    // Click handler
    overlayBtn.on('click', async function(event) {
        event.stopPropagation();

        if (isOverlayed) {
            await removeTrackOverlay(trackId);
        } else if (canAdd) {
            await addTrackOverlay(trackId);
        }
    });
}
```

### 6. Filter Reactivity Hook

**Location:** `gosling_mvp/frontend/d3_viewer.html` (in ClinVar filter checkbox handler ~line 2061)

```javascript
// In the ClinVar filter checkbox change handler:
cb.addEventListener('change', async (e) => {
    const labelKey = e.target.dataset.clinvarLabel;
    clinvarLabelFilter[labelKey] = e.target.checked;

    // Re-render track visualization
    if (selectedTracks.includes('clinvar.clinvar_label_list')) {
        renderChart();
    }

    // Update 3D spheres if ClinVar track is overlayed
    if (isTrackOverlayed('clinvar.clinvar_label_list')) {
        await updateTrackOverlaySpheres('clinvar.clinvar_label_list');
    }
});
```

### 7. Structure Change Handler

**Location:** `gosling_mvp/frontend/d3_viewer.html` (where structure is loaded)

```javascript
// When loading a new structure, clear existing overlays
async function loadStructure(gene) {
    // Clear existing track overlays
    await clearTrackOverlays();

    // ... existing structure loading code ...

    // Re-render chart to update overlay button states
    renderChart();
}
```

---

## File Changes Summary

### `gosling_mvp/frontend/src/molstar-plugin.ts`

**Add after line 34 (after existing constants):**
- `TRACK_OVERLAY_TAG_PREFIX` constant

**Add at end of file (after `disableDefaultBehaviors`):**
- `createTrackOverlaySpheres()` function
- `removeTrackOverlaySpheres()` function
- `clearAllTrackOverlays()` function
- Export all new functions

### `gosling_mvp/frontend/d3_viewer.html`

**Add after line ~1800 (after existing state variables):**
- `trackOverlayState` object
- `TRACK_OVERLAY_COLORS` constant
- `TRACK_OVERLAY_COLORS_CSS` constant
- `isEligibleForSphereOverlay()` function
- `isTrackOverlayed()` function
- `canAddTrackOverlay()` function
- `getNextColorSlot()` function
- `getTrackOverlayColor()` function
- `getTrackOverlayColorCSS()` function

**Add after state management functions:**
- `getResiduesWithTrackData()` function
- `addTrackOverlay()` function
- `removeTrackOverlay()` function
- `updateTrackOverlaySpheres()` function
- `clearTrackOverlays()` function

**Modify `renderChart()` function (~line 3023-3061):**
- Add overlay button after flip button for eligible tracks
- Adjust flip button position when overlay button is present

**Modify ClinVar filter handler (~line 2061):**
- Add call to `updateTrackOverlaySpheres()` when filter changes

**Modify structure loading function:**
- Add call to `clearTrackOverlays()` before loading new structure

---

## Testing Plan

### Unit Tests

1. **State Management**
   - `isEligibleForSphereOverlay()` returns true for ClinVar, Training, Constraint tracks
   - `isEligibleForSphereOverlay()` returns false for continuous tracks (O/E, counts)
   - `canAddTrackOverlay()` returns true when < 2 tracks overlayed
   - `getNextColorSlot()` returns correct slot based on current state

2. **Data Extraction**
   - `getResiduesWithTrackData()` returns correct residues for ClinVar data
   - ClinVar filter is respected in residue extraction
   - Positions without aa_pos are skipped
   - Positions from other genes are skipped

### Integration Tests

1. **Basic Overlay**
   - Load SCN2A structure
   - Add ClinVar track to selected tracks
   - Click overlay button
   - Verify blue spheres appear on structure

2. **Two Overlays**
   - Add first track overlay (blue)
   - Add second track overlay (cyan)
   - Verify both visible
   - Verify third track shows disabled button

3. **Remove Overlay**
   - Click overlay button on overlayed track
   - Verify spheres removed
   - Verify color slot freed for reuse

4. **Filter Reactivity**
   - Overlay ClinVar track
   - Toggle Pathogenic filter off
   - Verify pathogenic residue spheres removed
   - Toggle back on
   - Verify spheres restored

5. **Yellow Highlight Coexistence**
   - Overlay a track with spheres
   - Click on residue in 3D viewer
   - Verify yellow highlight appears on top of track spheres

6. **Structure Change**
   - Overlay tracks on one structure
   - Load different structure
   - Verify overlays cleared

### Manual Testing Checklist

- [ ] Button only appears for eligible tracks
- [ ] Button appears below flip button, both centered
- [ ] Button shows correct color when active
- [ ] Button disabled (grayed) when 2 tracks already overlayed
- [ ] Tooltip shows correct message for each state
- [ ] Spheres appear in correct positions
- [ ] Sphere size smaller than yellow halo
- [ ] Yellow halo visible when clicking overlayed residue
- [ ] Filter changes update spheres immediately
- [ ] No console errors during operations
- [ ] Performance acceptable with many spheres

---

## Verification Steps

1. **Start the application:**
   ```bash
   cd gosling_mvp/frontend && npm run dev
   cd gosling_mvp && python backend.py
   ```

2. **Open browser to http://localhost:5173**

3. **Load a gene with structure** (e.g., search for SCN2A)

4. **Add ClinVar track** to selected tracks (left panel tree)

5. **Verify overlay button** appears in track's colorbar area

6. **Click overlay button** - verify blue spheres on structure

7. **Add Constraint track** and overlay - verify cyan spheres

8. **Toggle ClinVar filter** - verify spheres update

9. **Click on structure** - verify yellow highlight coexists

10. **Remove overlay** - verify spheres disappear

11. **Build for production:**
    ```bash
    cd gosling_mvp/frontend && npm run build
    ```
