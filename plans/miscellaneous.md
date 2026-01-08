# Miscellaneous Feature Updates

## 1. Yellow Sphere Highlight - Click Outside to Remove

### Current Behavior
When a user clicks on a residue in the 3D structure viewer, a semi-transparent yellow sphere (spacefill) is rendered around the CA atom to highlight the selection. The sphere persists until another residue is clicked.

### Desired Behavior
Clicking outside of a residue (on empty space in the viewer, or outside the structure panel entirely) should remove the yellow sphere highlight and clear the selection.

### Implementation Location
- **File:** `browser/frontend/index.html`
- **Key functions:**
  - `createSphereHalo()` (line ~6410) - creates the sphere
  - `removeHalo()` (line ~6120) - removes the sphere
  - Click subscription (line ~6658) - handles residue clicks

### Implementation Approach

1. **Option A: Handle clicks on empty space within Mol* viewer**

   Modify the existing click subscription to detect when a click doesn't resolve to a valid residue:
   ```javascript
   plugin.behaviors.interaction.click.subscribe(async (e) => {
       const { current } = e;
       const loci = current?.loci;
       const seqId = getResidueFromLoci(loci);

       if (seqId !== null) {
           // Existing logic - create sphere
       } else {
           // NEW: Click on empty space - remove highlight
           await removeHalo();
           selectedResidue = null;
           lastClickedLoci = null;
           clearSelectionHighlight();
       }
   });
   ```

2. **Option B: Add document-level click listener for clicks outside structure panel**

   ```javascript
   document.addEventListener('click', (e) => {
       const viewerEl = document.getElementById('molstar-viewer');
       const structurePanel = document.getElementById('structure-panel');

       if (!viewerEl?.contains(e.target) && !structurePanel?.contains(e.target)) {
           removeHalo();
           selectedResidue = null;
           lastClickedLoci = null;
           clearSelectionHighlight();
       }
   });
   ```

3. **Recommended: Combine both approaches**
   - Option A handles clicks on empty space within the 3D viewer
   - Option B handles clicks completely outside the structure panel

### State Variables to Clear
```javascript
selectedResidue = null;      // Currently selected residue number
lastClickedLoci = null;      // Stored loci for re-applying halo
sphereHaloRef = null;        // Reference to sphere in state tree (cleared by removeHalo())
```

### Testing
- Click on residue → sphere appears
- Click on empty space in viewer → sphere disappears
- Click outside structure panel (e.g., on tracks panel) → sphere disappears
- Click on another residue → old sphere replaced with new one (existing behavior)

---

## 2. RGC Tracks Reorganization

### Current Problem
The RGC tracks are deeply nested, requiring too many clicks to access desired tracks. Current structure:

```
RGC
├── Summary
│   ├── Any
│   │   ├── Observed
│   │   ├── Expected μ
│   │   └── Max AF
│   ├── Synonymous (same structure)
│   └── Missense (same structure)
├── O/E Ratios
│   ├── Missense
│   │   ├── 3bp
│   │   │   ├── O/E
│   │   │   │   ├── (0)
│   │   │   │   ├── (10⁻⁶)
│   │   │   │   └── ...
│   │   │   └── Expected
│   │   │       └── ...
│   │   ├── 9bp (same structure)
│   │   ├── 21bp (same structure)
│   │   ├── 45bp (same structure)
│   │   └── 93bp (same structure)
│   ├── Synonymous (same structure)
│   └── Any (same structure)
└── VIRs
    ├── Mis
    │   ├── (0)
    │   │   ├── Length
    │   │   ├── Depth
    │   │   ├── Expected μ
    │   │   └── Mean Expected
    │   ├── (10⁻⁶) (same structure)
    │   └── ...
    ├── Syn (same structure)
    └── Any (same structure)
```

**Issues:**
- O/E Ratios: 5 levels deep (RGC → O/E Ratios → Consequence → Window → O/E|Expected → AF)
- VIRs: 4 levels deep (RGC → VIRs → Consequence → AF → Metric)
- Most users likely want a specific window size or AF threshold, not to browse all options

### Implementation Location
- **File:** `browser/backend/track_tree.py`
- **Key functions:**
  - `build_oe_tree()` (lines 10-67)
  - `build_vir_tree()` (lines 70-143)
  - Main tree assembly (lines 272-278)

### Proposed Reorganization

#### Option A: Flatten by Most Common Use Cases

Group by what users typically want to compare:

```
RGC
├── Counts
│   ├── Any
│   ├── Synonymous
│   └── Missense
├── Summary (existing)
├── O/E by Window (most used first)
│   ├── 21bp O/E
│   │   ├── Missense (0)
│   │   ├── Missense (10⁻⁴)
│   │   ├── Synonymous (0)
│   │   └── Any (0)
│   ├── 45bp O/E (same structure)
│   ├── 93bp O/E (same structure)
│   ├── 9bp O/E
│   └── 3bp O/E
├── O/E Expected (collapsed by default)
│   └── ... (same structure as O/E)
└── VIRs
    ├── Length
    │   ├── Missense (0)
    │   ├── Missense (10⁻⁴)
    │   ├── Synonymous (0)
    │   └── Any (0)
    ├── Depth (same structure)
    ├── Expected μ
    └── Mean Expected
```

**Benefits:**
- O/E tracks: 3 levels (RGC → Window → Consequence+AF)
- VIR tracks: 3 levels (RGC → Metric → Consequence+AF)
- Most common window sizes (21bp, 45bp, 93bp) listed first

#### Option B: Group by AF Threshold First

If users typically work at a specific AF threshold:

```
RGC
├── Counts
├── Summary
├── Strict (AF=0)
│   ├── O/E
│   │   ├── Missense [21bp, 45bp, 93bp]
│   │   ├── Synonymous [21bp, 45bp, 93bp]
│   │   └── Any [21bp, 45bp, 93bp]
│   └── VIRs
│       ├── Missense [Length, Depth, μ]
│       └── ...
├── Relaxed (AF≤10⁻⁴)
│   └── ... (same structure)
└── Very Relaxed (AF≤10⁻³)
    └── ...
```

#### Option C: Quick Access + Full Tree

Keep detailed tree but add quick-access shortcuts at the top:

```
RGC
├── ★ Quick Access
│   ├── O/E Mis 21bp (0)
│   ├── O/E Mis 45bp (0)
│   ├── O/E Syn 21bp (0)
│   ├── VIR Length Mis (0)
│   └── VIR Depth Mis (0)
├── Counts
├── Summary
├── O/E Ratios (full tree, collapsed by default)
└── VIRs (full tree, collapsed by default)
```

### Implementation Steps

1. **Modify `build_oe_tree()` in `track_tree.py`:**
   ```python
   def build_oe_tree() -> Dict[str, Any]:
       """Build flattened O/E tree grouped by window size."""
       windows = ["21bp", "45bp", "93bp", "9bp", "3bp"]  # Most used first
       consequences = [("mis", "Missense"), ("syn", "Synonymous"), ("any", "Any")]
       af_levels = [("0", "(0)"), ("1e-4", "(10⁻⁴)")]  # Most used AF levels

       return {
           "label": "O/E Ratios",
           "children": [
               {
                   "label": f"{window} O/E",
                   "children": [
                       {
                           "label": f"{cons_label} {af_label}",
                           "fieldId": f"rgc_{cons}_exomes_XX_XY_{window}_oe_{af_suffix}",
                           "percentileFieldId": f"dbnsfp.max_rgc_{cons}..._oe_{af_suffix}_exome_perc"
                       }
                       for cons, cons_label in consequences
                       for af_suffix, af_label in af_levels
                   ]
               }
               for window in windows
           ]
       }
   ```

2. **Modify `build_vir_tree()` similarly:**
   ```python
   def build_vir_tree() -> Dict[str, Any]:
       """Build flattened VIR tree grouped by metric."""
       metrics = ["Length", "Depth", "Expected μ", "Mean Expected"]
       # ... similar flattening logic
   ```

3. **Update frontend if needed** for any new grouping behavior

### Testing
- Verify all tracks still accessible
- Count clicks required to reach commonly-used tracks (target: ≤3 clicks)
- Ensure track field IDs unchanged (data binding preserved)
- Test track selection/deselection still works

---

## 3. 3D Panel Show/Hide Tab

### Current Behavior
The 3D structure panel has an "X" button to close it. Once closed, there's no intuitive way to reopen it.

### Desired Behavior
Replace the "X" button with a collapsible tab similar to the track selection panel. The tab should:
- Allow toggling the 3D panel open/closed
- Provide a visual indicator of the panel state
- Be consistent with other collapsible panels in the UI

### Implementation Location
- **File:** `browser/frontend/index.html`
- Look for structure panel container and close button

### Implementation Approach

1. **Remove the existing X button**

2. **Add a toggle tab** (similar to track panel):
   ```html
   <div id="structure-panel-toggle" class="panel-toggle">
     <span class="toggle-icon">◀</span>
     <span class="toggle-label">3D Structure</span>
   </div>
   ```

3. **CSS for toggle tab:**
   ```css
   .panel-toggle {
     position: absolute;
     top: 50%;
     transform: translateY(-50%);
     writing-mode: vertical-rl;
     text-orientation: mixed;
     background: #f0f0f0;
     padding: 10px 5px;
     cursor: pointer;
     border-radius: 0 4px 4px 0;
     z-index: 100;
   }

   .panel-toggle:hover {
     background: #e0e0e0;
   }

   .panel-toggle .toggle-icon {
     display: block;
     text-align: center;
     margin-bottom: 5px;
   }

   /* When panel is collapsed */
   .structure-panel.collapsed + .panel-toggle .toggle-icon {
     transform: rotate(180deg);
   }
   ```

4. **Toggle logic:**
   ```javascript
   document.getElementById('structure-panel-toggle').addEventListener('click', () => {
     const panel = document.getElementById('structure-panel');
     panel.classList.toggle('collapsed');

     // Update icon direction
     const icon = document.querySelector('#structure-panel-toggle .toggle-icon');
     icon.textContent = panel.classList.contains('collapsed') ? '▶' : '◀';

     // Trigger resize for other panels to fill space
     window.dispatchEvent(new Event('resize'));
   });
   ```

### Testing
- Click toggle tab → panel collapses, tab shows expand icon
- Click again → panel expands, tab shows collapse icon
- Panel state persists correctly during interactions
- Other panels resize appropriately when 3D panel toggles

---

## 4. Track Reordering on Track Panel

### Current Behavior
Track reordering is done in a separate region at the top of the track selection pane, away from where tracks are displayed.

### Desired Behavior
Allow drag-and-drop reordering directly on the track panel where the tracks are visualized. Users should be able to:
- Drag a track row up/down to reorder
- See visual feedback during drag
- Have changes reflected immediately in the visualization

### Implementation Location
- **File:** `browser/frontend/index.html`
- Track panel rendering and `selectedTracks` array management

### Implementation Approach

1. **Add drag handles to each track row:**
   ```html
   <div class="track-row" draggable="true" data-track-id="rgc_any_count">
     <span class="drag-handle">⋮⋮</span>
     <span class="track-label">Any Count</span>
     <span class="track-remove">×</span>
   </div>
   ```

2. **CSS for drag interaction:**
   ```css
   .track-row {
     display: flex;
     align-items: center;
     padding: 4px 8px;
     border-bottom: 1px solid #eee;
     cursor: default;
   }

   .drag-handle {
     cursor: grab;
     padding: 0 8px;
     color: #999;
   }

   .drag-handle:active {
     cursor: grabbing;
   }

   .track-row.dragging {
     opacity: 0.5;
     background: #f0f0f0;
   }

   .track-row.drag-over {
     border-top: 2px solid #007bff;
   }
   ```

3. **Drag-and-drop handlers:**
   ```javascript
   let draggedTrack = null;

   function initTrackDragAndDrop() {
     const trackRows = document.querySelectorAll('.track-row');

     trackRows.forEach(row => {
       row.addEventListener('dragstart', (e) => {
         draggedTrack = row;
         row.classList.add('dragging');
         e.dataTransfer.effectAllowed = 'move';
       });

       row.addEventListener('dragend', () => {
         row.classList.remove('dragging');
         document.querySelectorAll('.drag-over').forEach(el =>
           el.classList.remove('drag-over'));
         draggedTrack = null;
       });

       row.addEventListener('dragover', (e) => {
         e.preventDefault();
         if (row !== draggedTrack) {
           row.classList.add('drag-over');
         }
       });

       row.addEventListener('dragleave', () => {
         row.classList.remove('drag-over');
       });

       row.addEventListener('drop', (e) => {
         e.preventDefault();
         if (row !== draggedTrack) {
           // Reorder in selectedTracks array
           const fromId = draggedTrack.dataset.trackId;
           const toId = row.dataset.trackId;
           reorderTracks(fromId, toId);
           renderChart();
         }
       });
     });
   }

   function reorderTracks(fromId, toId) {
     const fromIdx = selectedTracks.indexOf(fromId);
     const toIdx = selectedTracks.indexOf(toId);

     if (fromIdx !== -1 && toIdx !== -1) {
       selectedTracks.splice(fromIdx, 1);
       selectedTracks.splice(toIdx, 0, fromId);
       updateSelectedTracksList();
     }
   }
   ```

4. **Remove old reordering UI** from track selection pane (or keep as alternative)

### Testing
- Drag track A above track B → order changes in visualization
- Drag track to top/bottom of list → works correctly
- Visual feedback during drag is clear
- Order persists after other interactions

---

## 5. Remove Title Bar and "Genes in View" Bar

### Current Behavior
The browser has:
- A title bar at the top
- A "Genes in view:" bar showing genes currently visible

### Desired Behavior
Remove both elements to maximize screen real estate for the actual data visualization.

### Implementation Location
- **File:** `browser/frontend/index.html`
- Look for title bar and genes-in-view elements

### Implementation Approach

1. **Identify elements to remove:**
   ```javascript
   // Find and comment out or remove:
   // - Title bar element (likely has id like "title-bar" or class "header")
   // - Genes in view bar (likely has id like "genes-in-view" or similar)
   ```

2. **Option A: Remove entirely**
   ```javascript
   // In HTML, remove or comment out:
   // <div id="title-bar">...</div>
   // <div id="genes-in-view">Genes in view: ...</div>
   ```

3. **Option B: Hide via CSS (easier to toggle back)**
   ```css
   #title-bar,
   #genes-in-view {
     display: none;
   }
   ```

4. **Adjust layout** to reclaim the vertical space:
   ```css
   /* Update container heights/positions to fill the space */
   #main-container {
     top: 0;  /* Was offset by title bar height */
   }
   ```

5. **Remove any JavaScript** that updates the "Genes in view" text to avoid unnecessary computation

### Considerations
- If title bar contains important controls (e.g., gene search), relocate them
- Consider making this a user preference/toggle
- Ensure removal doesn't break any layout calculations

### Testing
- Title bar and genes-in-view bar no longer visible
- Track panel and 3D panel use the reclaimed vertical space
- No JavaScript errors from removed elements
- Any relocated controls still functional

---

## Priority

1. **Yellow sphere click-outside** - Quick win, improves UX immediately
2. **Remove title/genes bars** - Quick win, more screen space
3. **3D panel show/hide tab** - Moderate effort, improves panel management
4. **Track reordering on panel** - Moderate effort, improves workflow
5. **RGC track reorganization** - Requires more design discussion on preferred structure

## Notes

- All changes are independent and can be implemented separately
- Track reorganization only affects the UI tree structure, not underlying data
- Consider adding user preferences for default expanded/collapsed state
- Features 3-5 are UI polish items that improve overall user experience
