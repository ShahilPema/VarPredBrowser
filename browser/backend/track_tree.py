"""
Track Tree Definitions for VarPredBrowser

Defines the hierarchical track structure for the genome browser frontend.
"""

from typing import Dict, Any, List


def build_oe_tree() -> Dict[str, Any]:
    """Build the O/E Ratios tree section."""

    def window_node(consequence: str, window: str, with_af: bool = False) -> Dict[str, Any]:
        if with_af:
            afs = [
                ("(0)", "af0epos00"),
                ("(10⁻⁶)", "af1eneg06"),
                ("(10⁻⁴)", "af1eneg04"),
            ]
        else:
            afs = [("(0)", "af0epos00")]
        return {
            "label": window,
            "children": [
                {
                    "label": "O/E",
                    "children": [
                        {
                            "label": lab,
                            "fieldId": f"rgc_{consequence}_exomes_XX_XY_{window}_oe_{suffix}",
                            "percentileFieldId": f"dbnsfp.max_rgc_{consequence}_exomes_XX_XY_{window}_oe_{suffix}_exome_perc"
                        }
                        for lab, suffix in afs
                    ],
                },
                {
                    "label": "Expected",
                    "children": [
                        {
                            "label": lab,
                            "fieldId": f"rgc_{consequence}_exomes_XX_XY_{window}_e_{suffix}",
                            "percentileFieldId": f"dbnsfp.max_rgc_{consequence}_exomes_XX_XY_{window}_e_{suffix}_exome_perc"
                        }
                        for lab, suffix in afs
                    ],
                },
            ],
        }

    miss_windows = [
        window_node("mis", "3bp", False),
        window_node("mis", "9bp", False),
        window_node("mis", "21bp", True),
        window_node("mis", "45bp", True),
        window_node("mis", "93bp", True),
    ]
    syn_windows = [window_node("syn", bp, False) for bp in ["3bp", "9bp", "21bp", "45bp", "93bp"]]
    any_windows = [window_node("any", bp, False) for bp in ["3bp", "9bp", "21bp", "45bp", "93bp"]]

    return {
        "label": "O/E Ratios",
        "children": [
            {"label": "Missense", "children": miss_windows},
            {"label": "Synonymous", "children": syn_windows},
            {"label": "Any", "children": any_windows},
        ],
    }


def build_vir_tree() -> Dict[str, Any]:
    """Build the VIRs tree section."""
    af_levels = [
        ("(0)", "af0epos00"),
        ("(10⁻⁶)", "af1eneg06"),
        ("(10⁻⁵)", "af1eneg05"),
        ("(10⁻⁴)", "af1eneg04"),
        ("(10⁻³)", "af1eneg03"),
        ("(10⁻²)", "af1eneg02"),
    ]

    def vir_group(consequence: str) -> Dict[str, Any]:
        return {
            "label": consequence.capitalize() if consequence != "any" else "Any",
            "children": [
                {
                    "label": lab,
                    "children": [
                        {
                            "label": "Length",
                            "fieldId": f"rgc_{consequence}_exomes_XX_XY_vir_length_{suffix}",
                            "percentileFieldId": f"dbnsfp.max_rgc_{consequence}_exomes_XX_XY_vir_length_{suffix}_exome_perc"
                        },
                        {
                            "label": "Depth",
                            "fieldId": f"rgc_{consequence}_exomes_XX_XY_vir_depth_{suffix}",
                            "percentileFieldId": f"dbnsfp.max_rgc_{consequence}_exomes_XX_XY_vir_depth_{suffix}_exome_perc"
                        },
                        {
                            "label": "Expected μ",
                            "fieldId": f"rgc_{consequence}_exomes_XX_XY_vir_mu_exp_{suffix}",
                            "percentileFieldId": f"dbnsfp.max_rgc_{consequence}_exomes_XX_XY_vir_mu_exp_{suffix}_exome_perc"
                        },
                        {
                            "label": "Mean Expected",
                            "fieldId": f"rgc_{consequence}_exomes_XX_XY_mean_vir_exp_{suffix}",
                            "percentileFieldId": f"dbnsfp.max_rgc_{consequence}_exomes_XX_XY_mean_vir_exp_{suffix}_exome_perc"
                        },
                    ],
                }
                for lab, suffix in af_levels
            ],
        }

    aa_children = [
        {
            "label": lab,
            "children": [
                {
                    "label": "Depth",
                    "fieldId": f"rgc_mis_exomes_XX_XY_aa_vir_depth_{suffix}",
                    "percentileFieldId": f"dbnsfp.max_rgc_mis_exomes_XX_XY_aa_vir_depth_{suffix}_exome_perc"
                },
                {
                    "label": "Length",
                    "fieldId": f"rgc_mis_exomes_XX_XY_aa_vir_length_{suffix}",
                    "percentileFieldId": f"dbnsfp.max_rgc_mis_exomes_XX_XY_aa_vir_length_{suffix}_exome_perc"
                },
            ],
        }
        for lab, suffix in af_levels
    ]

    miss = vir_group("mis")
    miss["children"] = miss["children"] + [{"label": "AA Level", "children": aa_children}]

    return {
        "label": "VIRs",
        "children": [
            miss,
            vir_group("syn"),
            vir_group("any"),
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
                                    {"label": "Max AF", "fieldId": "rgc_any_max_af"},
                                ],
                            },
                            {
                                "label": "Synonymous",
                                "children": [
                                    {"label": "Observed", "fieldId": "rgc_syn_obs_exomes_XX_XY"},
                                    {"label": "Expected μ", "fieldId": "rgc_syn_prob_mu_exomes_XX_XY"},
                                    {"label": "Max AF", "fieldId": "rgc_syn_max_af"},
                                ],
                            },
                            {
                                "label": "Missense",
                                "children": [
                                    {"label": "Observed", "fieldId": "rgc_mis_obs_exomes_XX_XY"},
                                    {"label": "Expected μ", "fieldId": "rgc_mis_prob_mu_exomes_XX_XY"},
                                    {"label": "Max AF", "fieldId": "rgc_mis_max_af"},
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
                    {"label": "Labels (Stacked)", "fieldId": "clinvar.clinvar_label_list", "type": "clinvar_stacked"},
                    {"label": "Variant Count", "fieldId": "clinvar.clinvar_count"},
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
                            {"label": "phyloP 17-way", "fieldId": "phylop17way"},
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
        'phylop17way': 'phyloP 17way',
        'ESM1b_score': 'ESM1b',
        'domain_name': 'Domains',
        'domain_id_interpro': 'Domain IDs',
        'domain_type': 'Domain Types',
        'source_db': 'Domain Sources',
        'domains': 'Protein Domains',
        # ClinVar columns
        'clinvar.clinvar_count': 'ClinVar Count',
        'clinvar.clinvar_label_list': 'ClinVar Labels',
        'clinvar.clinvar_status_list': 'ClinVar Status',
        'clinvar.clinvar_var_type_list': 'ClinVar Var Types',
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
    }
    if track_name in name_mappings:
        return name_mappings[track_name]

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
    if track_name.startswith('clinvar.'):
        return 'ClinVar'

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
    'mis_count_gt0': {
        'name': 'Missense Possible',
        'description': 'Positions where missense variants are possible (mis_count > 0)',
    },
    'any_count_gt0': {
        'name': 'Any Variant Possible',
        'description': 'All coding positions where any variant is possible (any_count > 0)',
    }
}

# Set of constraint stacked track fields that store variant arrays
CONSTRAINT_STACKED_FIELDS = {'Constraint', 'Core', 'Complete'}

# Set of dbNSFP stacked track fields that store variant arrays (allele, score, percentile)
DBNSFP_STACKED_FIELDS = {'AlphaMissense_stacked', 'ESM1b_stacked'}
