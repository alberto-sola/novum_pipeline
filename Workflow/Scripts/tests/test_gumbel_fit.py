from __future__ import annotations
import math
import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))  # Workflow/Scripts

import _gumbel_fit as gf

ALLOWED_IMPORTS = {"import math", "from collections import namedtuple"}


def test_module_is_dependency_free():
    # rule rnacalibrate runs in rnahybrid.yaml: no pandas, no numpy.
    source = Path(gf.__file__).read_text()
    imports = [l.strip() for l in source.splitlines()
               if l.startswith("import ") or l.startswith("from ")]
    extra = sorted(set(imports) - ALLOWED_IMPORTS)
    assert not extra, f"_gumbel_fit must stay stdlib-only, found: {extra}"


# --- Tier 1: hard rejects, straight from the RNAcalibrate source ---

def test_classify_fit_accepts_a_sound_fit():
    assert gf.classify_fit(233, 2.4792, 0.1673) is None

@pytest.mark.parametrize("sample_size", [-4, 0, 1])
def test_classify_fit_rejects_the_zeros_branch_by_its_source_condition(sample_size):
    # numerical.c:139 is `if (end<=start)`, i.e. N <= 1. A guard written as N < 0 from the
    # observed -4 would pass N=0 and N=1, which both emit 0.000000 0.000000.
    assert gf.classify_fit(sample_size, 0.0, 0.0) == "degenerate_sample"

def test_classify_fit_rejects_nan_before_comparing_it():
    # NaN <= 0 is False, so the sign checks alone would pass a NaN theta.
    assert gf.classify_fit(222, float("nan"), float("nan")) == "non_finite"
    assert gf.classify_fit(222, 2.4, float("inf")) == "non_finite"

def test_classify_fit_rejects_finite_but_physically_invalid_fits():
    # All three are measured outputs that pass an isnan guard untouched.
    assert gf.classify_fit(36, 13970.0, -2022.0) == "non_positive_theta"
    assert gf.classify_fit(10, -2070.0, 300.0) == "non_positive_xi"
    assert gf.classify_fit(44, -3.33, 0.803) == "non_positive_xi"

def test_classify_fit_rejects_non_finite_even_paired_with_a_negative_partner():
    # A non-finite field paired with a genuinely negative partner must still read as
    # non_finite; a reversed check order (sign checks before isfinite) would misclassify
    # these as non_positive_theta / non_positive_xi instead.
    assert gf.classify_fit(50, float("nan"), -5.0) == "non_finite"
    assert gf.classify_fit(50, -5.0, float("nan")) == "non_finite"


# --- Tier 3: the per-miRNA power law ---

QUERY_NT = 22

def _point(anchor, alpha, log_c, query_nt=QUERY_NT, sample_size=250):
    theta = math.exp(log_c - alpha * math.log(math.log(anchor * query_nt)))
    return gf.FitPoint(anchor, query_nt, theta, sample_size)

def test_se_of_log_theta_is_the_measured_coefficient_over_root_n():
    assert gf.SE_COEFFICIENT == 1.053
    assert gf.se_log_theta(250) == pytest.approx(1.053 / math.sqrt(250))

def test_fit_alpha_recovers_a_planted_exponent_exactly():
    points = [_point(a, 1.557, 0.5) for a in (76, 150, 300, 600, 900, 1500)]
    fit = gf.fit_alpha(points)
    assert fit.alpha == pytest.approx(1.557, abs=1e-9)
    assert fit.log_c == pytest.approx(0.5, abs=1e-9)

def test_fit_alpha_weights_by_sample_size_not_uniformly():
    # Five points at N=2000 sit exactly on the curve; one outlier at N=20 is 3x off it.
    # Correct N/SE^2 weighting suppresses the low-N outlier; unweighted or inverted
    # (SE^2/N) weighting instead lets it dominate and flips alpha negative.
    points = [_point(a, 1.557, 0.5, sample_size=2000) for a in (76, 150, 300, 600, 900)]
    outlier = _point(1500, 1.557, 0.5, sample_size=20)
    points.append(outlier._replace(theta=outlier.theta * 3.0))

    def alpha_with_weight(weight_of):
        sum_w = sum_x = sum_y = sum_xx = sum_xy = 0.0
        for p in points:
            x = math.log(math.log(p.target_length * p.query_length))
            y = math.log(p.theta)
            w = weight_of(p)
            sum_w += w; sum_x += w * x; sum_y += w * y
            sum_xx += w * x * x; sum_xy += w * x * y
        denom = sum_w * sum_xx - sum_x * sum_x
        return -(sum_w * sum_xy - sum_x * sum_y) / denom

    fit = gf.fit_alpha(points)
    assert fit.alpha == pytest.approx(1.5210490317517371, abs=1e-9)

    unweighted = alpha_with_weight(lambda p: 1.0)
    inverted = alpha_with_weight(lambda p: gf.SE_COEFFICIENT ** 2 / p.sample_size)
    assert abs(fit.alpha - unweighted) > 1.0
    assert abs(fit.alpha - inverted) > 1.0

def test_predict_theta_round_trips_the_planted_curve():
    point = _point(600, 1.557, 0.5)
    assert gf.predict_theta(1.557, 0.5, 600, QUERY_NT) == pytest.approx(point.theta)

def test_residual_ratio_is_one_on_the_curve():
    point = _point(600, 1.557, 0.5)
    assert gf.residual_ratio(point, 1.557, 0.5) == pytest.approx(1.0)

def test_fit_alpha_robust_drops_a_planted_outlier_and_recovers_the_true_alpha():
    points = [_point(a, 1.557, 0.5) for a in (76, 150, 300, 600, 900, 1500)]
    points[3] = points[3]._replace(theta=points[3].theta * 1.5)   # z = 4.8
    fit = gf.fit_alpha_robust(points)
    assert [p.target_length for p in fit.dropped] == [600]
    assert [p.target_length for p in fit.kept] == [76, 150, 300, 900, 1500]
    assert fit.alpha == pytest.approx(1.557, abs=1e-9)

def test_fit_alpha_robust_keeps_a_point_inside_three_sigma():
    points = [_point(a, 1.557, 0.5) for a in (76, 150, 300, 600, 900, 1500)]
    two_sigma = math.exp(2 * gf.se_log_theta(250))       # ~1.142 at N=250
    points[3] = points[3]._replace(theta=points[3].theta * two_sigma)
    fit = gf.fit_alpha_robust(points)
    assert fit.dropped == []

def test_fit_alpha_robust_rejects_a_curve_that_does_not_fall_with_length():
    # theta MUST decrease with target length; alpha <= 0 means the fit is describing
    # something other than the null and every extrapolation off it would be nonsense.
    points = [_point(a, -0.5, 0.5) for a in (76, 150, 300, 600, 900, 1500)]
    with pytest.raises(gf.UnusableCurveError, match="alpha"):
        gf.fit_alpha_robust(points)

def test_fit_alpha_robust_raises_when_fewer_than_three_points_are_given():
    # This trips fit_alpha's own upfront len(points) < MIN_POINTS guard, before the
    # rejection loop ever runs -- see the post-rejection variant below for that branch.
    with pytest.raises(gf.UnusableCurveError, match="at least 3"):
        gf.fit_alpha_robust([_point(76, 1.557, 0.5), _point(150, 1.557, 0.5)])

def test_fit_alpha_robust_raises_when_rejection_leaves_fewer_than_three_points():
    # 4 points at huge N (tight 3-sigma band); two are perturbed just enough to be
    # dropped, leaving 2 kept -- the post-rejection kept < MIN_POINTS branch.
    points = [_point(a, 1.557, 0.5, sample_size=100_000) for a in (76, 150, 300, 600)]
    points[2] = points[2]._replace(theta=points[2].theta * 1.02)
    points[3] = points[3]._replace(theta=points[3].theta * 0.98)
    with pytest.raises(gf.UnusableCurveError, match="survived"):
        gf.fit_alpha_robust(points)

def test_fit_alpha_reproduces_the_measured_third_mirna_ladder():
    # hsa-miR-1224-5p, 19 nt, k=10000 (spec 3.2). alpha was fitted on two OTHER miRNAs, so
    # this is out-of-sample: it must land near 1.557 with a ~5% residual spread.
    measured = [(76, 0.2532, 276), (150, 0.2285, 271), (300, 0.2024, 259),
                (600, 0.1742, 252), (900, 0.1673, 233), (1500, 0.1478, 245)]
    points = [gf.FitPoint(length, 19, theta, n) for length, theta, n in measured]
    fit = gf.fit_alpha_robust(points)
    assert fit.alpha == pytest.approx(1.5545, abs=0.01)
    assert fit.residual_spread == pytest.approx(0.0503, abs=0.005)
    assert fit.dropped == []          # nothing on this ladder is a 3-sigma outlier


# --- extrapolating xi ---

def test_weighted_median_picks_the_value_at_half_the_weight():
    assert gf.weighted_median([1.0, 2.0, 3.0], [1, 1, 1]) == 2.0

def test_weighted_median_is_dominated_by_the_heavy_point():
    assert gf.weighted_median([1.0, 2.0, 3.0], [1, 100, 1]) == 2.0
    assert gf.weighted_median([1.0, 2.0, 3.0], [100, 1, 1]) == 1.0

def test_weighted_median_ignores_input_order():
    assert gf.weighted_median([3.0, 1.0, 2.0], [1, 1, 1]) == 2.0

def test_weighted_median_rejects_an_empty_input():
    with pytest.raises(ValueError):
        gf.weighted_median([], [])
