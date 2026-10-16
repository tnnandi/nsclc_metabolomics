"""
Source-data corrections 
    "Combien GeoC& Metabo Epigenetic 08_26_2026 AM an AH without identifiers
     Responders vs Non-Responders .xlsx"
is a data refresh of metabolomics.new.xlsx: identical sheets, identical 142 x 691
shape, identical column names and order. Six cell-level defects are repaired here:

    C79   duplicate patient ID (row 1 and row 13 both read C69)
    C40   'Date of Initial Diagnosis' = 12/99/2018 (invalid day)
    C20   'BMI' = 30,09 (European decimal comma)
    C36   'Cancer YES or NO' = Control on a fully staged cancer record
    C55   'Date of death' / 'Last Day of Follow-up' transposed (two entries)

All six were confirmed with the data owner. C36 and C55 were adjudicated during the
2026-09-15 data-quality review; the other three are inferred from
metabolomics.new.xlsx.

Every repair matches on the exact bad value, so this module is a no-op once the
defects are fixed at source. Post-repair assertions hard-fail rather than let a
silently-coerced NaN propagate into the analysis.

Apply immediately after the raw load / header flatten, and BEFORE clean_data().
"""

from pathlib import Path

import numpy as np
import pandas as pd

try:
    from tabulate import tabulate
except ImportError:  # tabulate is in ai_env, but do not hard-fail on it
    tabulate = None

ID_COL = 'Double Deidentifier'

# --- The four known defects --------------------------------------------------
# kind='repair'        -> value is corrected, audited
# kind='informational' -> value is left alone, audited only
#
# 'row' is used only where the ID column is itself the broken field; everywhere
# else the row is located by patient ID so the fix survives any row reordering.
CORRECTIONS = [
    # Correction 1: Duplicate Patient ID resolution (Patient C79 vs C69)
    # Both row 1 and row 13 were labeled 'C69' in the spreadsheet, but represent
    # two distinct subjects (age 43 vs 61, completely different metabolite profiles).
    # Cross-referencing against metabolomics.new.xlsx confirmed row 1 is patient C79.
    # Repairing row 1's identifier from 'C69' to 'C79' prevents duplicate key collisions.
    dict(
        col=ID_COL, row=1, match_id=None,
        bad='C69', good='C79', kind='repair',
        reason=('Rows 1 and 13 both read C69 but are distinct subjects '
                '(age 43 vs 61, different metabolite vectors); '
                'metabolomics.new.xlsx has C79 at row 1'),
    ),

    # Correction 2: Malformed Date of Initial Diagnosis (Patient C40)
    # The diagnosis date contained an invalid day '12/99/2018' (data-entry typo for 2018-12-09).
    # This non-date string forced pandas to parse the entire column as an 'object' dtype
    # and caused C40 to be excluded from survival models due to an unparseable Time to Event.
    # Repaired to '2018-12-09' as confirmed in metabolomics.new.xlsx.
    dict(
        col='Date of Initial Diagnosis', row=None, match_id='C40',
        bad='12/99/2018', good='2018-12-09', kind='repair',
        reason=('Invalid day 99 poisons the column to object dtype and drops '
                'this patient for missing time-to-event; '
                'metabolomics.new.xlsx has 2018-12-09'),
    ),

    # Correction 3: European Decimal Comma in BMI (Patient C20)
    # BMI was entered as '30,09' with a comma instead of a period. This prevented numeric
    # casting and poisoned the entire BMI column into an 'object' dtype.
    # Repaired to '30.09' to restore numeric parsing.
    dict(
        col='BMI', row=None, match_id='C20',
        bad='30,09', good='30.09', kind='repair',
        reason=('European decimal comma poisons the column to object dtype; '
                'metabolomics.new.xlsx has 30.09'),
    ),

    # Correction 4: Misclassified Cohort Label (Patient C36)
    # Patient C36 was flagged as 'Control' in 'Cancer YES or NO', despite possessing
    # a full cancer clinical profile (Squamous cell carcinoma, Stage I, T1 N1 M0, 0.8 cm tumor,
    # diagnosis 2019-10-04, surgery 2019-12-06). C36 was the sole non-cancer entry with tumor staging.
    # Confirmed as CANCER with the data owner; reclassified from 'Control' to 'CANCER' so the
    # patient is correctly included in the cancer survival cohort (alive / censored at 132 d).
    dict(
        col='Cancer YES or NO', row=None, match_id='C36',
        bad='Control', good='CANCER', kind='repair',
        reason=('Flagged Control in both sheets while carrying a complete cancer '
                'record: Squamous histology, Stage I, T1 N1 M0, 0.8 cm tumor, '
                'diagnosis 2019-10-04, surgery 2019-12-06. It is the only '
                'non-cancer row in the workbook with any staging, and '
                'metabolomics.new.xlsx has CANCER. Restored to CANCER per the data '
                'owner. C36 is Squamous so it clears the SCLC filter, and is alive, '
                'so it joins the survival cohort as a censored observation '
                '(TTE 132 d from baseline draw).'),
    ),

    # Correction 5: Transposed Date of Death (Patient C55 - Part 1)
    # Patient C55 was the only decedent whose recorded Date of Death (2021-11-11) occurred
    # BEFORE the Last Day of Follow-up (2021-12-04), representing an impossible chronological reversal.
    # Adjudicated with the data owner: the two date fields were accidentally transposed during entry.
    # Swapped Date of Death to '2021-12-04' (true date of death).
    dict(
        col='Date of death', row=None, match_id='C55',
        bad='2021-11-11 00:00:00', good='2021-12-04', kind='repair',
        reason=('C55 is the only deceased patient whose Date of death (2021-11-11) '
                'precedes Last Day of Follow-up (2021-12-04); the two dates are '
                'transposed. Both sheets carry the same pair, so this is one '
                'upstream entry error, not a transcription divergence. Swapped per '
                'the data owner. Nothing else on the row falls between the two '
                'dates (IO ended 2021-08-11), so the swap is internally consistent.'),
    ),

    # Correction 6: Transposed Last Follow-up Date (Patient C55 - Part 2)
    # Companion fix to Correction 5. Swapped Last Day of Follow-up to '2021-11-11'
    # so that the final clinical follow-up visit logically precedes the verified death date.
    dict(
        col='Last Day of Follow-up', row=None, match_id='C55',
        bad='2021-12-04 00:00:00', good='2021-11-11', kind='repair',
        reason='Second half of the C55 death / last-follow-up transposition (see above).',
    ),
]

# Columns the residual-dirt scan expects to be fully numeric after correction.
NUMERIC_CLINICAL = [
    'Age at the time of diagnosis', 'BMI', 'Pack Year History',
    'Size 1 of tumor (cms)', 'Age at time of surgery',
    'Vital Status\n0=Alive\n1=Dead',
]
METHYLATION_COLS = ['CDO1. P', 'TAC1. P', 'SOX17.P', 'HOXA7.', 'HOXA9.P', 'ZFP42.P']
SOCIAL_COLS = [
    'ADI', 'EKI', 'PMC', 'RHI', 'TBC', 'CCG', 'CRV', 'EDB', 'EKG', 'EKR', 'FAI',
    'HBS', 'HDX', 'HTA', 'ICR', 'INC', 'LNG', 'MCD', 'NVC', 'POV', 'UMP', 'UNS',
]
DATE_COLS = ['Date of Initial Diagnosis', 'Date of death', 'Last Day of Follow-up']

# Strings that legitimately mean "missing" in this workbook.
NA_STRINGS = {'', 'NA', 'N/A', 'NAN', 'NONE', '-', 'NULL', '#N/A', 'UNKNOWN', 'UNK'}

# Metabolites are selected as Baseline_C* minus these substrings; the four
# Baseline_C00669_* follow-up timing flags are Yes/No and must stay excluded.
_MET_EXCLUDE = ['Date', 'Follow Up', 'Draw', 'Immunotherapy', 'Chemotherapy', 'SBRT']


def _is_missing(v):
    return (v is None or v is pd.NaT or (isinstance(v, float) and np.isnan(v))
            or (isinstance(v, str) and v.strip().upper() in NA_STRINGS))


def _analysis_metabolites(df):
    return [c for c in df.columns
            if c.startswith('Baseline_C') and not any(k in c for k in _MET_EXCLUDE)]


def apply_corrections(df, audit_path=None, verbose=True):
    """Repair the known defects, audit every change, and assert the result is clean.

    Returns (df, audit_df). Raises AssertionError if a correction failed to take
    or if an unexpected defect remains.
    """
    df = df.copy()
    audit = []

    if verbose:
        print("\n--- Applying source-data corrections ---")

    for c in CORRECTIONS:
        col, bad, good, kind = c['col'], c['bad'], c['good'], c['kind']

        if col not in df.columns:
            raise AssertionError(f"correction target column missing: {col!r}")

        # Locate the row: positionally when the ID itself is broken, else by ID.
        if c['row'] is not None:
            idx = [df.index[c['row']]]
        else:
            ids = df[ID_COL].astype(str).str.strip()
            idx = list(df.index[ids == c['match_id']])

        # Guard the C69 repair: only act if the duplicate pre-state is what we saw.
        if col == ID_COL and kind == 'repair':
            ids = df[ID_COL].astype(str).str.strip()
            dups = sorted(df.index[ids == bad].tolist())
            if len(dups) < 2:
                if verbose:
                    print(f"  [skip] {ID_COL}: {bad!r} is not duplicated "
                          f"(already fixed upstream)")
                continue
            if dups != [df.index[1], df.index[13]]:
                raise AssertionError(
                    f"unexpected {bad!r} duplicate rows {dups}; expected rows 1 and 13. "
                    "Source layout changed - re-audit before trusting this repair.")

        for i in idx:
            cur = df.at[i, col]
            cur_s = str(cur).strip()
            pid = str(df.at[i, ID_COL]).strip()

            if cur_s != bad:
                if verbose:
                    print(f"  [skip] {col!r} @ {pid}: found {cur_s!r}, "
                          f"expected {bad!r} (already fixed upstream)")
                continue

            if kind == 'repair':
                value = good
                if pd.api.types.is_datetime64_any_dtype(df[col]):
                    # A str assignment here would demote the column to object
                    # dtype and turn every NaT into the string 'NaT'.
                    value = pd.Timestamp(good)
                df.at[i, col] = value
                new_pid = good if col == ID_COL else pid
                if verbose:
                    print(f"  [fix ] {col!r} @ {pid}: {bad!r} -> {good!r}")
            else:
                new_pid = pid
                if verbose:
                    print(f"  [note] {col!r} @ {pid}: {bad!r} left as-is")

            audit.append(dict(patient_id=new_pid, column=col, old_value=bad,
                              new_value=good if kind == 'repair' else bad,
                              kind=kind, reason=c['reason']))

    # --- Hard assertions: the whole point of the module ----------------------
    ids = df[ID_COL].astype(str).str.strip()
    dupes = [k for k, v in ids.value_counts().items() if v > 1]
    assert not dupes, f"duplicate patient IDs after correction: {dupes}"

    bmi = pd.to_numeric(df['BMI'], errors='coerce')
    bad_bmi = [(ids[i], repr(v)) for i, v in df['BMI'].items()
               if pd.isna(bmi[i]) and not _is_missing(v)]
    assert not bad_bmi, f"BMI still non-numeric after correction: {bad_bmi}"

    for dc in DATE_COLS:
        if dc not in df.columns:
            continue
        parsed = pd.to_datetime(df[dc], errors='coerce')
        bad_dt = [(ids[i], repr(v)) for i, v in df[dc].items()
                  if pd.isna(parsed[i]) and not _is_missing(v)]
        assert not bad_dt, f"{dc!r} still unparseable after correction: {bad_dt}"

    for group, cols in (('methylation', METHYLATION_COLS),
                        ('social', SOCIAL_COLS),
                        ('metabolite', _analysis_metabolites(df))):
        present = [c for c in cols if c in df.columns]
        dirty = {c: [repr(v) for v in df[c] if isinstance(v, str) and not _is_missing(v)][:3]
                 for c in present}
        dirty = {k: v for k, v in dirty.items() if v}
        assert not dirty, f"non-numeric values in {group} columns: {dirty}"

    # Death must never precede the last documented follow-up: the follow-up window
    # has to extend to the terminal event. Survival durations take the death date
    # for Event=1, so a transposed pair does not corrupt the duration itself, but
    # it does signal a bad row -- fail loudly rather than model it.
    if {'Date of death', 'Last Day of Follow-up'} <= set(df.columns):
        death = pd.to_datetime(df['Date of death'], errors='coerce')
        lastfu = pd.to_datetime(df['Last Day of Follow-up'], errors='coerce')
        inverted = [(ids[i], str(death[i].date()), str(lastfu[i].date()))
                    for i in df.index
                    if pd.notna(death[i]) and pd.notna(lastfu[i]) and death[i] < lastfu[i]]
        assert not inverted, (
            "Date of death precedes Last Day of Follow-up for "
            f"(id, death, last_followup): {inverted}")

    if verbose:
        print(f"  [ok  ] {len(ids)} unique patient IDs; BMI, dates, "
              f"methylation, social and metabolite columns all clean; "
              f"no death-before-last-follow-up inversions")

    # --- Residual-dirt scan: surface anything new, loudly --------------------
    _residual_scan(df, ids, verbose=verbose)

    audit_df = pd.DataFrame(audit, columns=['patient_id', 'column', 'old_value',
                                            'new_value', 'kind', 'reason'])
    if audit_path is not None:
        audit_path = Path(audit_path)
        audit_path.parent.mkdir(parents=True, exist_ok=True)
        audit_df.to_csv(audit_path, index=False)
        if verbose:
            print(f"  Audit trail -> {audit_path}")

    if verbose and len(audit_df):
        cols = ['patient_id', 'column', 'old_value', 'new_value', 'kind']
        if tabulate is not None:
            print(tabulate(audit_df[cols], headers='keys', tablefmt='github',
                           showindex=False))
        else:
            print(audit_df[cols].to_string(index=False))
        print()

    return df, audit_df


def _residual_scan(df, ids, verbose=True):
    """Report any non-numeric / unparseable value not covered by CORRECTIONS."""
    import re

    known = {(c['match_id'] or str(c['bad'])) + '|' + c['col'] for c in CORRECTIONS}
    findings = []

    for col in NUMERIC_CLINICAL:
        if col not in df.columns:
            continue
        for i, v in df[col].items():
            if isinstance(v, str) and not _is_missing(v) \
               and not re.fullmatch(r'-?\d+(\.\d+)?', v.strip()):
                if f"{ids[i]}|{col}" not in known:
                    findings.append((ids[i], col, repr(v), 'non-numeric'))

    for col in DATE_COLS:
        if col not in df.columns:
            continue
        for i, v in df[col].items():
            if isinstance(v, str) and not _is_missing(v) \
               and pd.isna(pd.to_datetime(v, errors='coerce')):
                if f"{ids[i]}|{col}" not in known:
                    findings.append((ids[i], col, repr(v), 'unparseable date'))

    if findings:
        print("\n  !! RESIDUAL DATA DEFECTS not covered by CORRECTIONS:")
        for pid, col, val, kind in findings:
            print(f"     {pid}  {col!r}  {val}  ({kind})")
        print("     Add these to CORRECTIONS or confirm them with the data providers.\n")
    elif verbose:
        print("  [ok  ] residual-dirt scan found no uncatalogued defects")

    return findings


def provenance_banner(file_path, extra=None):
    """Print input provenance so any figure traces to the exact input revision."""
    import datetime as _dt

    p = Path(file_path)
    print("=" * 78)
    print("INPUT PROVENANCE")
    print(f"  file  : {p.name}")
    print(f"  path  : {p.resolve()}")
    if p.exists():
        st = p.stat()
        mtime = _dt.datetime.fromtimestamp(st.st_mtime).strftime('%Y-%m-%d %H:%M:%S')
        print(f"  mtime : {mtime}   size: {st.st_size:,} bytes")
    if extra:
        for k, v in extra.items():
            print(f"  {k:6}: {v}")
    print("=" * 78)
