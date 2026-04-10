#!/usr/bin/env python3

import fnmatch


REQUIRED_KEYS = [
    'name',
    'treatments',
    'controls',
    'good_barcodes_file',
]

ALLOWED_KEYS = {
    'name',
    'treatments',
    'controls',
    'good_barcodes_file',
    'treatments_negative_selection',
    'controls_negative_selection',
}


def as_string_list(value):
    if value is None:
        return []
    if isinstance(value, list):
        return [str(v) for v in value]
    return [str(value)]


def match_selectors(selectors, all_ids):
    matched = set()
    explicit_not_found = []

    for selector in selectors:
        if '*' in selector or '?' in selector:
            matched.update(fnmatch.filter(all_ids, selector))
        else:
            if selector in all_ids:
                matched.add(selector)
            else:
                explicit_not_found.append(selector)

    return sorted(matched), sorted(explicit_not_found)


def resolve_selectors(selectors, all_ids, negative_selectors=None):
    negative_selectors = negative_selectors or []

    positive_ids, explicit_not_found = match_selectors(selectors, all_ids)
    negative_ids, negative_explicit_not_found = match_selectors(negative_selectors, all_ids)

    resolved = sorted(set(positive_ids) - set(negative_ids))
    excluded_by_negative_selection = sorted(set(positive_ids) & set(negative_ids))

    return {
        'resolved_ids': resolved,
        'explicit_not_found': explicit_not_found,
        'negative_explicit_not_found': negative_explicit_not_found,
        'excluded_by_negative_selection': excluded_by_negative_selection,
    }


def derive_status(treat_ids, control_ids):
    if treat_ids and control_ids:
        return 'OK'
    if not treat_ids and not control_ids:
        return 'ALL_SAMPLES_MISSING'
    return 'SOME_SAMPLES_MISSING'


def validate_comparison_schema(comparison):
    keys = set(comparison.keys())
    invalid = sorted(keys - ALLOWED_KEYS)
    missing = [k for k in REQUIRED_KEYS if k not in keys]

    if missing:
        raise ValueError(
            "Missing required fields in comparison '{}': {}. Required fields are: {}. You probably have a typo in your comparisons.json :)".format(
                comparison.get('name', '<unnamed>'),
                ', '.join(missing),
                ', '.join(REQUIRED_KEYS),
            )
        )

    if invalid:
        raise ValueError(
            "Invalid fields in comparison '{}': {}. Allowed fields are: {}. You probably have a typo in your comparisons.json :)".format(
                comparison.get('name', '<unnamed>'),
                ', '.join(invalid),
                ', '.join(sorted(ALLOWED_KEYS)),
            )
        )


def build_comparison_resolution(comparison, all_ids):
    validate_comparison_schema(comparison)

    name = str(comparison['name'])
    treat_sel = as_string_list(comparison.get('treatments'))
    ctrl_sel = as_string_list(comparison.get('controls'))
    treat_neg_sel = as_string_list(comparison.get('treatments_negative_selection'))
    ctrl_neg_sel = as_string_list(comparison.get('controls_negative_selection'))

    treat_info = resolve_selectors(treat_sel, all_ids, treat_neg_sel)
    ctrl_info = resolve_selectors(ctrl_sel, all_ids, ctrl_neg_sel)

    treat_ids = treat_info['resolved_ids']
    control_ids = ctrl_info['resolved_ids']

    treat_ids_set = set(treat_ids)
    control_ids_set = set(control_ids)
    treat_neg_set = set(treat_neg_sel)
    control_neg_set = set(ctrl_neg_sel)

    # Preserve selector order in status detail while matching legacy semantics.
    missing_treat = [sid for sid in treat_sel if sid not in treat_ids_set and sid not in treat_neg_set]
    missing_control = [sid for sid in ctrl_sel if sid not in control_ids_set and sid not in control_neg_set]

    if not missing_treat and not missing_control:
        status = 'OK'
        status_detail = ''
    elif (treat_sel and not treat_ids) or (ctrl_sel and not control_ids):
        status = 'ALL_SAMPLES_MISSING'
        parts = []
        if missing_treat:
            parts.append('treatments: ' + ', '.join(missing_treat))
        if missing_control:
            parts.append('controls: ' + ', '.join(missing_control))
        status_detail = ' | '.join(parts)
    else:
        status = 'SOME_SAMPLES_MISSING'
        parts = []
        if missing_treat:
            parts.append('treatments: ' + ', '.join(missing_treat))
        if missing_control:
            parts.append('controls: ' + ', '.join(missing_control))
        status_detail = ' | '.join(parts)

    return {
        'name': name,
        'good_barcodes_file': comparison['good_barcodes_file'],
        'treatment_selectors': treat_sel,
        'control_selectors': ctrl_sel,
        'treatment_negative_selectors': treat_neg_sel,
        'control_negative_selectors': ctrl_neg_sel,
        'treat_info': treat_info,
        'control_info': ctrl_info,
        'treat_ids': treat_ids,
        'control_ids': control_ids,
        'status': status,
        'status_detail': status_detail,
    }


def build_all_comparison_resolutions(comparisons_payload, all_ids):
    comparisons = comparisons_payload.get('comparisons', [])
    return [build_comparison_resolution(cmp, all_ids) for cmp in comparisons]
