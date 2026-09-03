"""
Mondrian (class-conditional) split conformal prediction for binary cytotoxicity classification.

This file intentionally contains ONLY the conformal-prediction logic.
It does not split data and does not train any model.

Expected workflow in the model-building code
--------------------------------------------
1. Create TRAIN / CALIBRATION / TEST splits.
2. Fit a classifier using TRAIN only.
3. Obtain P(y=1) on CALIBRATION and TEST.
4. Pass those probabilities to `run_conformal_binary` below.

Label convention used here:
    0 = non-cytotoxic
    1 = cytotoxic

Nonconformity score:
    s(x, c) = 1 - P(c | x)

A candidate class c is included in the conformal prediction set when:
    s(x, c) <= q_c

where q_c is the class-specific threshold estimated from calibration samples
whose TRUE label is c.
"""

import numpy as np
import pandas as pd


# -----------------------------------------------------------------------------
# 1. Build class-specific conformal thresholds from the calibration set
# -----------------------------------------------------------------------------

def mondrian_thresholds_binary(y_cal, p_cyto_cal, alpha_by_class=None):
    """
    Learn class-specific Mondrian conformal thresholds.

    Parameters
    ----------
    y_cal : array-like, shape (n_cal,)
        True calibration labels. Must contain 0/1 labels.

    p_cyto_cal : array-like, shape (n_cal,)
        Model-predicted P(y=1) = P(cytotoxic) on the calibration set.

    alpha_by_class : dict, optional
        Class-specific error tolerances.
        Default: {0: 0.10, 1: 0.10}, corresponding to nominal 90% coverage
        for both classes.

        Example for greater protection of the cytotoxic class:
            {0: 0.10, 1: 0.05}
        which targets 90% coverage for class 0 and 95% for class 1.

    Returns
    -------
    score_thresholds : dict
        Nonconformity thresholds q_c for classes 0 and 1.

    threshold_table : pandas.DataFrame
        Human-readable threshold details, including equivalent probability
        cutoffs P(c) >= 1 - q_c.
    """
    if alpha_by_class is None:
        alpha_by_class = {0: 0.10, 1: 0.10}

    y_cal = np.asarray(y_cal, dtype=int)
    p_cyto_cal = np.asarray(p_cyto_cal, dtype=float)

    if y_cal.ndim != 1 or p_cyto_cal.ndim != 1:
        raise ValueError("y_cal and p_cyto_cal must be one-dimensional arrays.")
    if len(y_cal) != len(p_cyto_cal):
        raise ValueError("y_cal and p_cyto_cal must have the same length.")
    if not np.all(np.isin(y_cal, [0, 1])):
        raise ValueError("y_cal must contain only binary labels 0 and 1.")
    if np.any(~np.isfinite(p_cyto_cal)) or np.any((p_cyto_cal < 0) | (p_cyto_cal > 1)):
        raise ValueError("Predicted probabilities must be finite values in [0, 1].")

    # Columns are P(class 0) and P(class 1).
    probs_cal = np.column_stack([1.0 - p_cyto_cal, p_cyto_cal])

    score_thresholds = {}
    rows = []

    for cls in (0, 1):
        mask = y_cal == cls
        n_cls = int(mask.sum())
        if n_cls == 0:
            raise ValueError(f"No calibration examples were found for class {cls}.")

        alpha = float(alpha_by_class[cls])
        if not 0.0 < alpha < 1.0:
            raise ValueError(f"alpha for class {cls} must lie strictly between 0 and 1.")

        # For each calibration example of TRUE class c:
        # nonconformity = 1 - probability assigned to its true class.
        scores = 1.0 - probs_cal[mask, cls]
        scores = np.sort(scores)

        # Finite-sample split-conformal rank (1-indexed).
        k = int(np.ceil((n_cls + 1) * (1.0 - alpha)))

        # If the requested coverage is too high for the available number of
        # calibration samples, the conservative conformal threshold is +inf.
        if k > n_cls:
            q = np.inf
            probability_cutoff = 0.0
        else:
            q = float(scores[k - 1])
            probability_cutoff = float(1.0 - q)

        score_thresholds[cls] = q

        rows.append({
            "class": cls,
            "n_calibration": n_cls,
            "alpha": alpha,
            "target_coverage": 1.0 - alpha,
            "conformal_rank": k,
            "score_threshold": q,
            "equivalent_probability_cutoff": probability_cutoff,
        })

    return score_thresholds, pd.DataFrame(rows)


# -----------------------------------------------------------------------------
# 2. Apply the thresholds to test/new samples
# -----------------------------------------------------------------------------

def predict_conformal_sets_binary(p_cyto, score_thresholds):
    """
    Construct conformal prediction sets for test/new samples.

    Parameters
    ----------
    p_cyto : array-like, shape (n_samples,)
        Model-predicted P(y=1) = P(cytotoxic).

    score_thresholds : dict
        Thresholds returned by `mondrian_thresholds_binary`.

    Returns
    -------
    prediction_sets : list[set]
        Possible outputs for binary classification:
            {0}     -> non-cytotoxic only
            {1}     -> cytotoxic only
            {0, 1}  -> both classes plausible / uncertain
            set()   -> neither class included
    """
    p_cyto = np.asarray(p_cyto, dtype=float)

    if p_cyto.ndim != 1:
        raise ValueError("p_cyto must be a one-dimensional array.")
    if np.any(~np.isfinite(p_cyto)) or np.any((p_cyto < 0) | (p_cyto > 1)):
        raise ValueError("Predicted probabilities must be finite values in [0, 1].")

    probs = np.column_stack([1.0 - p_cyto, p_cyto])
    scores = 1.0 - probs

    prediction_sets = []

    for i in range(len(p_cyto)):
        current_set = set()

        for cls in (0, 1):
            if scores[i, cls] <= score_thresholds[cls]:
                current_set.add(cls)

        prediction_sets.append(current_set)

    return prediction_sets


# -----------------------------------------------------------------------------
# 3. Evaluate conformal prediction on an untouched labelled test set
# -----------------------------------------------------------------------------

def evaluate_conformal_binary(y_test, prediction_sets, alpha_by_class=None):
    """
    Evaluate coverage, efficiency and cytotoxicity-specific safety metrics.

    Parameters
    ----------
    y_test : array-like
        True test labels.

    prediction_sets : list[set]
        Conformal sets returned by `predict_conformal_sets_binary`.

    alpha_by_class : dict, optional
        Used only to report the target class-wise coverage.

    Returns
    -------
    metrics : dict
    """
    if alpha_by_class is None:
        alpha_by_class = {0: 0.10, 1: 0.10}

    y_test = np.asarray(y_test, dtype=int)

    if len(y_test) != len(prediction_sets):
        raise ValueError("y_test and prediction_sets must have the same length.")
    if not np.all(np.isin(y_test, [0, 1])):
        raise ValueError("y_test must contain only binary labels 0 and 1.")

    covered = np.array([
        true_label in pred_set
        for true_label, pred_set in zip(y_test, prediction_sets)
    ])

    set_sizes = np.array([len(pred_set) for pred_set in prediction_sets])

    metrics = {
        "coverage": float(covered.mean()),
        "singleton_rate": float(np.mean(set_sizes == 1)),
        "doubleton_rate": float(np.mean(set_sizes == 2)),
        "empty_rate": float(np.mean(set_sizes == 0)),
        "avg_set_size": float(set_sizes.mean()),
    }

    # Class-specific coverage.
    for cls in (0, 1):
        mask = y_test == cls
        metrics[f"coverage_class_{cls}"] = (
            float(covered[mask].mean()) if mask.any() else np.nan
        )
        metrics[f"target_coverage_class_{cls}"] = 1.0 - float(alpha_by_class[cls])

    # Safety metrics, assuming class 1 = cytotoxic.
    cyto_indices = np.where(y_test == 1)[0]

    if len(cyto_indices) > 0:
        cyto_sets = [prediction_sets[i] for i in cyto_indices]

        # True cytotoxic sample was labelled ONLY non-cytotoxic.
        metrics["false_safe_rate"] = float(
            np.mean([pred_set == {0} for pred_set in cyto_sets])
        )

        # Cytotoxic class is absent from the conformal set entirely.
        metrics["cytotoxic_exclusion_rate"] = float(
            np.mean([1 not in pred_set for pred_set in cyto_sets])
        )
    else:
        metrics["false_safe_rate"] = np.nan
        metrics["cytotoxic_exclusion_rate"] = np.nan

    return metrics


# -----------------------------------------------------------------------------
# 4. Convenience wrapper: calibration -> thresholds -> test sets -> evaluation
# -----------------------------------------------------------------------------

def run_conformal_binary(
    y_cal,
    p_cyto_cal,
    y_test,
    p_cyto_test,
    alpha_by_class=None,
):
    """
    Run the complete conformal-prediction stage for one already-trained model.

    Returns
    -------
    threshold_table : pandas.DataFrame
    prediction_sets : list[set]
    metrics : dict
    """
    if alpha_by_class is None:
        alpha_by_class = {0: 0.10, 1: 0.10}

    score_thresholds, threshold_table = mondrian_thresholds_binary(
        y_cal=y_cal,
        p_cyto_cal=p_cyto_cal,
        alpha_by_class=alpha_by_class,
    )

    prediction_sets = predict_conformal_sets_binary(
        p_cyto=p_cyto_test,
        score_thresholds=score_thresholds,
    )

    metrics = evaluate_conformal_binary(
        y_test=y_test,
        prediction_sets=prediction_sets,
        alpha_by_class=alpha_by_class,
    )

    return threshold_table, prediction_sets, metrics


# -----------------------------------------------------------------------------
# EXAMPLE OF HOW TO CALL THIS FROM THE MODEL-BUILDING CODE
# -----------------------------------------------------------------------------
#
# For a sklearn-style model:
#
#     model.fit(X_train, y_train)
#
#     p_cal  = model.predict_proba(X_cal)[:, 1]   # P(cytotoxic)
#     p_test = model.predict_proba(X_test)[:, 1]  # P(cytotoxic)
#
#     alpha = {0: 0.10, 1: 0.10}
#
#     thresholds, conformal_sets, metrics = run_conformal_binary(
#         y_cal=y_cal,
#         p_cyto_cal=p_cal,
#         y_test=y_test,
#         p_cyto_test=p_test,
#         alpha_by_class=alpha,
#     )
#
#     print(thresholds)
#     print(pd.DataFrame([metrics]))
#
# To protect the cytotoxic class more strongly later:
#
#     alpha = {0: 0.10, 1: 0.05}
#
# IMPORTANT:
# - TRAIN is used to fit the model.
# - CALIBRATION is used to estimate conformal thresholds.
# - TEST is used only for final evaluation.
# - Do not use calibration samples to train the model whose probabilities are
#   used for conformal calibration.
