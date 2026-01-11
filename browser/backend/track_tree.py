"""
Track Tree Definitions for VarPredBrowser

Defines the hierarchical track structure for the genome browser frontend.
"""

from typing import Dict, Any, List


def build_oe_tree() -> Dict[str, Any]:
    """Build the O/E Ratios tree section - flattened by window size with Raw/%ile separation."""
    # Window sizes in ascending order
    windows = ["3bp", "9bp", "21bp", "45bp", "93bp"]
    consequences = [("mis", "Missense"), ("syn", "Synonymous"), ("any", "Any")]
    af_suffix = "af0epos00"  # Only AF≥0

    def build_window_section(window: str) -> Dict[str, Any]:
        raw_children = []
        perc_children = []
        for cons, cons_label in consequences:
            raw_children.append({
                "label": cons_label,
                "fieldId": f"rgc_{cons}_exomes_XX_XY_{window}_oe_{af_suffix}",
            })
            perc_children.append({
                "label": cons_label,
                "fieldId": f"rgc_{cons}_exomes_XX_XY_{window}_oe_{af_suffix}_exome_perc",
            })
        return {
            "label": f"{window} O/E",
            "children": [
                {"label": "Raw", "children": raw_children},
                {"label": "Exome-Wide %ile", "children": perc_children},
            ]
        }

    return {
        "label": "O/E Ratios",
        "children": [build_window_section(w) for w in windows],
    }


def build_vir_tree() -> Dict[str, Any]:
    """Build the VIRs tree section - flattened by metric type with Raw/%ile separation."""
    # AF levels with ≥ notation
    af_levels = [
        ("AF≥0", "af0epos00"),
        ("AF≥10⁻⁴", "af1eneg04"),
        ("AF≥10⁻⁶", "af1eneg06"),
    ]
    consequences = [("mis", "Missense"), ("syn", "Synonymous"), ("any", "Any")]
    # Metrics: (label, suffix, has_percentile)
    metrics = [
        ("Length", "vir_length", True),
        ("Depth", "vir_depth", False),  # No percentiles for depth
        ("Expected μ", "vir_mu_exp", True),
        ("Mean Expected", "mean_vir_exp", True),
    ]

    def build_metric_section(metric_label: str, metric_suffix: str, has_percentile: bool) -> Dict[str, Any]:
        raw_children = []
        perc_children = []
        for cons, cons_label in consequences:
            for af_label, af_suffix in af_levels:
                raw_children.append({
                    "label": f"{cons_label} {af_label}",
                    "fieldId": f"rgc_{cons}_exomes_XX_XY_{metric_suffix}_{af_suffix}",
                })
                if has_percentile:
                    perc_children.append({
                        "label": f"{cons_label} {af_label}",
                        "fieldId": f"rgc_{cons}_exomes_XX_XY_{metric_suffix}_{af_suffix}_exome_perc",
                    })
        children = [{"label": "Raw", "children": raw_children}]
        if has_percentile:
            children.append({"label": "Exome-Wide %ile", "children": perc_children})
        return {"label": metric_label, "children": children}

    # Build AA Level section (missense only) - depth has no percentiles, length does
    aa_raw_depth = []
    aa_raw_length = []
    aa_perc_length = []
    for af_label, af_suffix in af_levels:
        aa_raw_depth.append({
            "label": f"Depth {af_label}",
            "fieldId": f"rgc_mis_exomes_XX_XY_aa_vir_depth_{af_suffix}",
        })
        aa_raw_length.append({
            "label": f"Length {af_label}",
            "fieldId": f"rgc_mis_exomes_XX_XY_aa_vir_length_{af_suffix}",
        })
        aa_perc_length.append({
            "label": f"Length {af_label}",
            "fieldId": f"rgc_mis_exomes_XX_XY_aa_vir_length_{af_suffix}_exome_perc",
        })

    aa_section = {
        "label": "AA Level (Missense)",
        "children": [
            {"label": "Raw", "children": aa_raw_depth + aa_raw_length},
            {"label": "Exome-Wide %ile", "children": aa_perc_length},
        ]
    }

    return {
        "label": "VIRs",
        "children": [
            build_metric_section(label, suffix, has_perc) for label, suffix, has_perc in metrics
        ] + [aa_section],
    }


def build_gnomad_oe_tree() -> Dict[str, Any]:
    """Build the gnomAD O/E Ratios tree section (Section 6 of browser-data-refactor.md)."""

    def window_node(consequence: str, window: str) -> Dict[str, Any]:
        return {
            "label": window,
            "children": [
                {
                    "label": "O/E",
                    "fieldId": f"gnomad_{consequence}_{window}_oe",
                },
                {
                    "label": "Expected",
                    "fieldId": f"gnomad_{consequence}_{window}_e",
                },
            ],
        }

    windows = ["3bp", "9bp", "21bp", "45bp", "93bp"]  # Ascending order
    consequences = [("mis", "Missense"), ("syn", "Synonymous"), ("any", "Any")]

    return {
        "label": "O/E Ratios",
        "children": [
            {
                "label": cons_label,
                "children": [window_node(cons, w) for w in windows],
            }
            for cons, cons_label in consequences
        ],
    }


def build_coverage_tree() -> Dict[str, Any]:
    """Build the Coverage track tree section (browser-data-refactor.md Section 2)."""
    return {
        "label": "Coverage",
        "children": [
            {
                "label": "Exomes",
                "children": [
                    {"label": "Mean", "fieldId": "gnomad_exomes_mean"},
                    {"label": "Median", "fieldId": "gnomad_exomes_median"},
                    {"label": "Over 20x", "fieldId": "gnomad_exomes_over_20"},
                    {"label": "Over 30x", "fieldId": "gnomad_exomes_over_30"},
                ],
            },
            {
                "label": "Genomes",
                "children": [
                    {"label": "Mean", "fieldId": "gnomad_genomes_mean"},
                    {"label": "Median", "fieldId": "gnomad_genomes_median"},
                    {"label": "Over 20x", "fieldId": "gnomad_genomes_over_20"},
                    {"label": "Over 30x", "fieldId": "gnomad_genomes_over_30"},
                ],
            },
        ],
    }


def build_variant_frequency_tree() -> Dict[str, Any]:
    """Build the Variant Frequencies track tree section."""
    return {
        "label": "Variant Frequencies",
        "children": [
            {"label": "gnomAD Exomes", "fieldId": "gnomad_exomes_variants", "type": "variant_frequency"},
            {"label": "gnomAD Genomes", "fieldId": "gnomad_genomes_variants", "type": "variant_frequency"},
            {"label": "RGC", "fieldId": "rgc_variants", "type": "variant_frequency"},
        ],
    }


def build_track_tree() -> Dict[str, Any]:
    """Build the complete track tree structure."""
    track_tree = {
        "label": "Tracks",
        "children": [
            # 1) Counts
            {
                "label": "Counts",
                "children": [
                    {"label": "Any", "fieldId": "rgc_any_count"},
                    {"label": "Synonymous", "fieldId": "rgc_syn_count"},
                    {"label": "Missense", "fieldId": "rgc_mis_count"},
                ],
            },
            # 2) RGC group
            {
                "label": "RGC",
                "children": [
                    {
                        "label": "Summary",
                        "children": [
                            {
                                "label": "Any",
                                "children": [
                                    {"label": "Observed", "fieldId": "rgc_any_obs_exomes_XX_XY"},
                                    {"label": "Expected μ", "fieldId": "rgc_any_prob_mu_exomes_XX_XY"},
                                ],
                            },
                            {
                                "label": "Synonymous",
                                "children": [
                                    {"label": "Observed", "fieldId": "rgc_syn_obs_exomes_XX_XY"},
                                    {"label": "Expected μ", "fieldId": "rgc_syn_prob_mu_exomes_XX_XY"},
                                ],
                            },
                            {
                                "label": "Missense",
                                "children": [
                                    {"label": "Observed", "fieldId": "rgc_mis_obs_exomes_XX_XY"},
                                    {"label": "Expected μ", "fieldId": "rgc_mis_prob_mu_exomes_XX_XY"},
                                ],
                            },
                        ],
                    },
                ],
            },
            # 3) ClinVar
            {
                "label": "ClinVar",
                "children": [
                    {"label": "Variants", "fieldId": "clinvar_variants", "type": "clinvar_variants"},
                ],
            },
            # 4) Training Labels
            {
                "label": "Training Labels",
                "children": [
                    {"label": "Labelled", "fieldId": "training.train_counts.labelled"},
                    {"label": "Unlabelled", "fieldId": "training.train_counts.unlabelled"},
                    {"label": "Labelled (High Qual)", "fieldId": "training.train_counts.labelled_high_qual"},
                    {"label": "Unlabelled (High Qual)", "fieldId": "training.train_counts.unlabelled_high_qual"},
                ],
            },
            # 5) Domains
            {
                "label": "Domains",
                "children": [
                    {"label": "All Domains", "fieldId": "domains", "type": "domain_array"},
                ],
            },
            # 6) Constraint Predictions
            {
                "label": "Constraint Predictions",
                "children": [
                    {"label": "Constraint_1000", "fieldId": "Constraint", "type": "constraint_stacked"},
                    {"label": "Core_1000", "fieldId": "Core", "type": "constraint_stacked"},
                    {"label": "Complete_1000", "fieldId": "Complete", "type": "constraint_stacked"},
                ],
            },
            # 7) Comparators
            {
                "label": "Comparators",
                "children": [
                    {
                        "label": "Conservation",
                        "children": [
                            {"label": "phyloP 447-way", "fieldId": "phylop_scores_447way"},
                            {"label": "phyloP 100-way", "fieldId": "phylop_scores_100way"},
                        ],
                    },
                    {
                        "label": "Constraint",
                        "children": [
                            {"label": "MTR (dbNSFP)", "fieldId": "dbnsfp.max_RGC_MTR_MTR"},
                            {"label": "MTR Percentile", "fieldId": "dbnsfp.max_RGC_MTR_MTRpercentile_exome"},
                            {"label": "CCR Residual Percentile", "fieldId": "dbnsfp.max_Non_Neuro_CCR_resid_pctile"},
                        ],
                    },
                    {
                        "label": "Pathogenicity / Model Scores",
                        "children": [
                            {"label": "AlphaMissense (max)", "fieldId": "dbnsfp.max_AlphaMissense_am_pathogenicity"},
                            {"label": "AlphaMissense (stacked)", "fieldId": "AlphaMissense_stacked", "type": "dbnsfp_stacked"},
                            {"label": "ESM1b (max)", "fieldId": "dbnsfp.max_ESM1b_score"},
                            {"label": "ESM1b (stacked)", "fieldId": "ESM1b_stacked", "type": "dbnsfp_stacked"},
                        ],
                    },
                    {
                        "label": "AlphaSync",
                        "children": [
                            {"label": "pLDDT", "fieldId": "dbnsfp.max_AlphaSync_plddt"},
                            {"label": "pLDDT (10-window)", "fieldId": "dbnsfp.max_AlphaSync_plddt10"},
                            {"label": "RelASA", "fieldId": "dbnsfp.max_AlphaSync_relasa"},
                            {"label": "RelASA (10-window)", "fieldId": "dbnsfp.max_AlphaSync_relasa10"},
                        ],
                    },
                ],
            },
        ],
    }

    # Add O/E and VIRs under RGC
    for child in track_tree["children"]:
        if child["label"] == "RGC":
            child["children"].append(build_oe_tree())
            child["children"].append(build_vir_tree())
            break

    # Add gnomAD section (Section 6 of browser-data-refactor.md)
    gnomad_section = {
        "label": "gnomAD",
        "children": [
            build_gnomad_oe_tree(),
        ],
    }
    # Insert gnomAD after RGC
    for i, child in enumerate(track_tree["children"]):
        if child["label"] == "RGC":
            track_tree["children"].insert(i + 1, gnomad_section)
            break

    # Add Coverage as top-level section after gnomAD
    for i, child in enumerate(track_tree["children"]):
        if child["label"] == "gnomAD":
            track_tree["children"].insert(i + 1, build_coverage_tree())
            break

    # Add Variant Frequencies section after Coverage
    for i, child in enumerate(track_tree["children"]):
        if child["label"] == "Coverage":
            track_tree["children"].insert(i + 1, build_variant_frequency_tree())
            break

    return track_tree


def simplify_track_name(track_name: str) -> str:
    """Simplify track name by removing common prefixes/terms."""
    # Direct mappings for new annotation tracks
    name_mappings = {
        'MTR_RGC': 'MTR',
        'Non_Neuro_CCR_resid_pctile': 'CCR',
        'AlphaMissense_am_pathogenicity': 'AlphaMissense',
        'phylop_scores_447way': 'phyloP 447way',
        'phylop_scores_100way': 'phyloP 100way',
        'ESM1b_score': 'ESM1b',
        'domain_name': 'Domains',
        'domain_id_interpro': 'Domain IDs',
        'domain_type': 'Domain Types',
        'source_db': 'Domain Sources',
        'domains': 'Protein Domains',
        # ClinVar columns
        'clinvar_variants': 'ClinVar Variants',
        'clinvar.clinvar_count': 'ClinVar Count',
        'clinvar.clinvar_label_list': 'ClinVar Labels',
        'clinvar.clinvar_status_list': 'ClinVar Status',
        'clinvar.clinvar_var_type_list': 'ClinVar Var Types',
        # Variant frequency tracks
        'rgc_variants': 'RGC Variants',
        'gnomad_exomes_variants': 'gnomAD Exome Variants',
        'gnomad_genomes_variants': 'gnomAD Genome Variants',
        # Training columns
        'training.train_counts.labelled': 'Labelled',
        'training.train_counts.unlabelled': 'Unlabelled',
        'training.train_counts.labelled_high_qual': 'Labelled (HQ)',
        'training.train_counts.unlabelled_high_qual': 'Unlabelled (HQ)',
        # dbNSFP columns
        'dbnsfp.max_AlphaMissense_am_pathogenicity': 'AlphaMissense',
        'dbnsfp.max_RGC_MTR_MTR': 'MTR',
        'dbnsfp.max_RGC_MTR_MTRpercentile_exome': 'MTR Percentile',
        'dbnsfp.max_Non_Neuro_CCR_resid_pctile': 'CCR',
        'dbnsfp.max_ESM1b_score': 'ESM1b',
        'dbnsfp.max_AlphaSync_plddt': 'pLDDT',
        'dbnsfp.max_AlphaSync_plddt10': 'pLDDT (10)',
        'dbnsfp.max_AlphaSync_relasa': 'RelASA',
        'dbnsfp.max_AlphaSync_relasa10': 'RelASA (10)',
        # Stacked tracks
        'AlphaMissense_stacked': 'AlphaMissense (stacked)',
        'ESM1b_stacked': 'ESM1b (stacked)',
        # gnomAD coverage
        'gnomad_exomes_mean': 'Exome Mean Cov',
        'gnomad_exomes_median': 'Exome Median Cov',
        'gnomad_exomes_over_20': 'Exome >20x',
        'gnomad_exomes_over_30': 'Exome >30x',
        'gnomad_genomes_mean': 'Genome Mean Cov',
        'gnomad_genomes_median': 'Genome Median Cov',
        'gnomad_genomes_over_20': 'Genome >20x',
        'gnomad_genomes_over_30': 'Genome >30x',
    }
    if track_name in name_mappings:
        return name_mappings[track_name]

    # Handle gnomAD constraint O/E columns
    if track_name.startswith('gnomad_') and ('_oe' in track_name or track_name.endswith('_e')):
        parts = track_name.replace('gnomad_', '').split('_')
        return f"gnomAD {' '.join(parts).upper()}"

    # Handle cross-norm percentile columns
    if '_cross_norm_perc' in track_name:
        base = track_name.replace('_cross_norm_perc', '')
        return f"{base.replace('_', ' ').title()} (cross-norm %)"

    simplified = track_name

    # Remove prefixes
    simplified = simplified.replace('dbnsfp.max_', '')
    simplified = simplified.replace('clinvar.', '')
    simplified = simplified.replace('training.train_counts.', '')
    simplified = simplified.replace('rgc_', '')
    simplified = simplified.replace('_exomes_', '_')
    simplified = simplified.replace('exomes_', '')
    simplified = simplified.replace('_XX_XY_', '_')
    simplified = simplified.replace('XX_XY_', '')
    simplified = simplified.replace('_af0epos00', '')
    simplified = simplified.replace('af0epos00', '')

    # Handle percentile suffix
    if simplified.endswith('_exome_perc'):
        simplified = simplified.replace('_exome_perc', ' (Percentile)')

    # Make it more readable
    simplified = simplified.replace('_', ' ').title()
    return simplified


def categorize_track(track_name: str) -> str:
    """Categorize track by its name into hierarchical categories."""
    name_lower = track_name.lower()

    # ClinVar annotations
    if track_name.startswith('clinvar.') or track_name == 'clinvar_variants':
        return 'ClinVar'

    # Variant frequency tracks
    elif track_name in ('rgc_variants', 'gnomad_exomes_variants', 'gnomad_genomes_variants'):
        return 'Variant Frequencies'

    # Training labels
    elif track_name.startswith('training.'):
        return 'Training Labels'

    # dbNSFP scores
    elif track_name.startswith('dbnsfp.'):
        if 'alphasync' in name_lower:
            return 'AlphaSync'
        elif 'alphamissense' in name_lower:
            return 'Pathogenicity'
        elif 'esm1b' in name_lower:
            return 'Pathogenicity'
        elif 'mtr' in name_lower:
            return 'Constraint'
        elif 'ccr' in name_lower:
            return 'Constraint'
        elif '_exome_perc' in name_lower:
            return 'Percentiles'
        else:
            return 'dbNSFP Scores'

    # Domain annotations
    elif 'domain' in name_lower or track_name == 'domains':
        return 'Domains'

    # Conservation scores (PhyloP)
    elif 'phylop' in name_lower:
        return 'Conservation'

    # gnomAD constraint O/E metrics
    elif track_name.startswith('gnomad_') and ('_oe' in name_lower or name_lower.endswith('_e')):
        return 'gnomAD Constraint'

    # gnomAD coverage
    elif track_name.startswith('gnomad_') and ('mean' in name_lower or 'median' in name_lower or 'over_' in name_lower):
        return 'gnomAD Coverage'

    # Cross-normalized percentiles
    elif '_cross_norm_perc' in name_lower:
        return 'Cross-Norm Percentiles'

    # dbNSFP stacked tracks
    elif track_name.endswith('_stacked'):
        return 'Pathogenicity (Stacked)'

    # Pathogenicity scores (AlphaMissense, ESM1b)
    elif 'alphamissense' in name_lower or 'esm1b' in name_lower:
        return 'Pathogenicity'

    # Constraint scores (MTR, CCR)
    elif 'mtr' in name_lower or 'ccr' in name_lower:
        return 'Constraint'

    # O/E Ratios by window size
    elif '_oe_' in name_lower:
        if '3bp' in name_lower:
            return 'O/E Ratios > 3bp Window'
        elif '93bp' in name_lower:
            return 'O/E Ratios > 93bp Window'
        elif '1000bp' in name_lower or '1kb' in name_lower:
            return 'O/E Ratios > 1kb Window'
        else:
            return 'O/E Ratios > Other Windows'

    # VIR metrics with subcategories
    elif 'vir' in name_lower:
        if 'length' in name_lower:
            return 'VIRs > VIR Length'
        elif 'depth' in name_lower:
            return 'VIRs > VIR Depth'
        elif 'mean' in name_lower and 'exp' in name_lower:
            return 'VIRs > Mean VIR Expression'
        elif 'mu' in name_lower and 'exp' in name_lower:
            return 'VIRs > Mu Expression'
        else:
            return 'VIRs > Other Metrics'

    # Counts
    elif 'count' in name_lower:
        return 'Counts'

    # Observations
    elif '_o_' in name_lower or 'obs' in name_lower:
        return 'Observations'

    # Expectations
    elif '_e_' in name_lower or 'exp' in name_lower:
        return 'Expectations'

    # Allele frequencies
    elif 'af' in name_lower or 'max_af' in name_lower:
        return 'Allele Frequencies'

    else:
        return 'RGC Data'


# Build the track tree once at module load
TRACK_TREE = build_track_tree()

# Available filters
FILTERS = {
    'missense_only': {
        'name': 'Missense Only',
        'description': 'Positions where missense variants are possible (mis_count > 0)',
    },
    'synonymous_only': {
        'name': 'Synonymous Only',
        'description': 'Positions where synonymous variants are possible (syn_count > 0)',
    },
    'all_sites': {
        'name': 'All Sites',
        'description': 'All coding positions where any variant is possible (any_count > 0)',
    }
}

# Set of constraint stacked track fields that store variant arrays
CONSTRAINT_STACKED_FIELDS = {'Constraint', 'Core', 'Complete'}

# Set of dbNSFP stacked track fields that store variant arrays (allele, score, percentile)
DBNSFP_STACKED_FIELDS = {'AlphaMissense_stacked', 'ESM1b_stacked'}

# Set of variant frequency track fields
VARIANT_FREQUENCY_FIELDS = {'rgc_variants', 'gnomad_exomes_variants', 'gnomad_genomes_variants'}

# ClinVar variants field (new struct format)
CLINVAR_VARIANTS_FIELD = 'clinvar_variants'
