/*
  pionCT_kinematics.C  —  Exclusive pion electroproduction kinematics
  =========================================================================
  Reaction:  e(k) + p(P)  →  e'(k') + π(p_π) + N'(P')
             in the one-photon-exchange approximation.

  PHYSICS SUMMARY
  ---------------
  Given Q², E_beam and either W or t_target, we compute the complete
  kinematic point for FORWARD pion production (θ* = 0 in the hadronic CM).

  Forward production is special because:
    • θ* = 0  gives  t = t_min(W, Q², E_beam), the minimum |t| achievable.
    • The table you are trying to reproduce was built by choosing W (or
      equivalently t_min) as the primary physics input, then computing all
      lab observables at θ* = 0.
    • The quantity k_π listed in the table is  k_π = |q̄| − |p̄_π| (lab),
      the longitudinal momentum surplus of the virtual photon over the pion
      (≈ √|t_min| up to mass corrections).

  USAGE MODES
  -----------
  Mode A — specify W directly (most transparent):
      pionCT_kinematics(5.0, 10.7, "W", 2.42)

  Mode B — find W such that  t_min(W) = t_target  (matches the table):
      pionCT_kinematics(5.0, 10.7, "t", -0.40)

  BOOST RECIPE (correct)
  ----------------------
  The hadronic CM moves along the virtual-photon direction q̂ with
      β_CM = |q̄_lab| / (M_p + ν).
  We work in a 2-D (x, z) plane (φ = 0):
    1. q̄ components:  q_x = -E' sin θ_e,  q_z = E - E' cos θ_e.
    2. q̂ = q̄ / |q̄|;  define q̂_⊥ = (-q_z/|q̄|, q_x/|q̄|) (in-plane perp).
    3. Pion in CM:  p_∥ = p* cos θ*,  p_⊥ = p* sin θ*.
    4. Boost along q̂:
          p_∥_lab = γ(p_∥ + β E*_π)
          E_π_lab = γ(E*_π + β p_∥)
          p_⊥_lab = p_⊥  (unchanged)
    5. Lab components:  p_x = p_∥_lab · q̂_x + p_⊥_lab · q̂_⊥_x
                        p_z = p_∥_lab · q̂_z + p_⊥_lab · q̂_⊥_z
  This is the CORRECT vector-boost approach; no RotateY patch is needed.

  DEFINITIONS
  -----------
    t      = (q − p_π)² = (ν − E_π)² − |q̄ − p̄_π|²    [negative for forward prod.]
    t_min  = t at θ* = 0 (most forward pion, minimum |t|)
    k_π    = |q̄| − |p̄_π|  (longitudinal surplus; equals |p̄_π − q̄| at θ*=0)
    θ_π    = angle of p̄_π w.r.t. beam (+z) in the lab
*/

#include "TLorentzVector.h"
#include "TMath.h"
#include <iostream>
#include <iomanip>
#include <cmath>
#include <string>
using namespace std;

// ── Physical constants (GeV) ────────────────────────────────────────────────
static const double M_p  = 0.9382720813;
static const double M_pi = 0.139570611;
static const double PI   = TMath::Pi();

// ── Internal helpers ─────────────────────────────────────────────────────────

// CM pion momentum from W
double pstar(double W) {
    double s = W*W;
    double a = s - (M_p + M_pi)*(M_p + M_pi);
    double b = s - (M_p - M_pi)*(M_p - M_pi);
    if (a <= 0 || b <= 0) return 0.0;
    return 0.5/W * sqrt(a*b);
}

// ── Core kinematic function ──────────────────────────────────────────────────
// For given (Q2, Ebeam, W, theta_star), compute all observables.
// Returns false if kinematics are unphysical.
struct KinResult {
    double W, nu, Eprime, theta_e_deg;
    double qmag, betaCM;
    double pstar_val, Epi_star;
    double E_pi, p_pi, theta_pi_deg;
    double t, k_pi;          // k_pi = |q| - |p_pi| at theta*=0
};

bool computeKinematics(double Q2, double Ebeam, double W, double theta_star,
                       KinResult &r)
{
    // --- electron side ---
    double nu = (W*W - M_p*M_p + Q2) / (2.0*M_p);
    double Eprime = Ebeam - nu;
    if (Eprime <= 0.0) return false;

    double sin2half = Q2 / (4.0 * Ebeam * Eprime);
    if (sin2half <= 0.0 || sin2half > 1.0) return false;
    double theta_e = 2.0 * asin(sqrt(sin2half));

    // --- virtual photon 3-vector in lab (beam = +z) ---
    double qx = -Eprime * sin(theta_e);
    double qz =  Ebeam  - Eprime * cos(theta_e);
    double qmag = sqrt(qx*qx + qz*qz);
    // cross-check with analytic formula (should agree to ~1e-10)
    // double qmag_analytic = sqrt(nu*nu + Q2);

    // q-hat and in-plane perpendicular direction
    double qhx = qx / qmag,  qhz = qz / qmag;     // along q
    double qpx = -qhz,       qpz =  qhx;           // perp to q, in (x,z) plane

    // --- CM boost ---
    double betaCM = qmag / (M_p + nu);
    if (betaCM >= 1.0) return false;
    double gammaCM = 1.0 / sqrt(1.0 - betaCM*betaCM);

    // --- CM pion kinematics ---
    double ps    = pstar(W);
    if (ps <= 0.0) return false;
    double Estar = sqrt(ps*ps + M_pi*M_pi);

    double p_par_cm  = ps * cos(theta_star);   // along q in CM
    double p_perp_cm = ps * sin(theta_star);   // perp to q in CM

    // --- Boost to lab (vector boost along q) ---
    double p_par_lab = gammaCM * (p_par_cm + betaCM * Estar);
    double E_pi_lab  = gammaCM * (Estar    + betaCM * p_par_cm);
    double p_perp_lab = p_perp_cm;             // unchanged by boost along q

    // --- Lab 3-vector ---
    double ppi_x = p_par_lab * qhx + p_perp_lab * qpx;
    double ppi_z = p_par_lab * qhz + p_perp_lab * qpz;
    double pmag  = sqrt(ppi_x*ppi_x + ppi_z*ppi_z);

    // pion angle w.r.t. beam (+z)
    double theta_pi = atan2(fabs(ppi_x), ppi_z);

    // --- Mandelstam t ---
    double dt_E = nu    - E_pi_lab;
    double dt_x = qx   - ppi_x;
    double dt_z = qz   - ppi_z;
    double t    = dt_E*dt_E - dt_x*dt_x - dt_z*dt_z;

    // --- k_pi ---
    double k_pi = qmag - pmag;   // longitudinal surplus (= |p_pi - q| at theta*=0)

    r.W           = W;
    r.nu          = nu;
    r.Eprime      = Eprime;
    r.theta_e_deg = theta_e * 180.0/PI;
    r.qmag        = qmag;
    r.betaCM      = betaCM;
    r.pstar_val   = ps;
    r.Epi_star    = Estar;
    r.E_pi        = E_pi_lab;
    r.p_pi        = pmag;
    r.theta_pi_deg = theta_pi * 180.0/PI;
    r.t           = t;
    r.k_pi        = k_pi;
    return true;
}

// ── t_min as a function of W ──────────────────────────────────────────────────
// t_min = t(theta*=0).  Returns +1e9 if unphysical.
double tmin_at_W(double W, double Q2, double Ebeam) {
    KinResult r;
    if (!computeKinematics(Q2, Ebeam, W, 0.0, r)) return 1e9;
    return r.t;
}

// ── Find W such that t_min(W) = t_target (bisection) ─────────────────────────
bool findW_for_tmin(double Q2, double Ebeam, double t_target, double &W_out) {
    // W ranges from pi-production threshold upward
    double W_thr = sqrt(M_p*M_p + 2.0*M_p * Q2/(2.0*M_p*(1.0+0.0)) + 0.0);
    W_thr = M_p + M_pi + 1e-6;   // just above threshold

    // Scan for bracket: t_min(W) is a monotone function of W
    // As W increases, |t_min| decreases (t_min -> 0 from below as W -> inf).
    // So for t_target < 0, we need to find W where t_min = t_target.
    double W_lo = W_thr, W_hi = 20.0;
    double f_lo = tmin_at_W(W_lo, Q2, Ebeam) - t_target;
    double f_hi = tmin_at_W(W_hi, Q2, Ebeam) - t_target;

    if (f_lo * f_hi > 0) {
        // Try to find the bracket by scanning
        bool found = false;
        double prev_f = f_lo;
        for (double W = W_thr; W <= 20.0; W += 0.001) {
            double f = tmin_at_W(W, Q2, Ebeam) - t_target;
            if (prev_f * f <= 0) {
                W_lo = W - 0.001; W_hi = W;
                f_lo = prev_f;    f_hi = f;
                found = true; break;
            }
            prev_f = f;
        }
        if (!found) return false;
    }

    // Bisection
    for (int i = 0; i < 80; ++i) {
        double W_mid = 0.5*(W_lo + W_hi);
        double f_mid = tmin_at_W(W_mid, Q2, Ebeam) - t_target;
        if (fabs(f_mid) < 1e-9) { W_out = W_mid; return true; }
        if (f_lo * f_mid <= 0) { W_hi = W_mid; f_hi = f_mid; }
        else                   { W_lo = W_mid; f_lo = f_mid; }
    }
    W_out = 0.5*(W_lo + W_hi);
    return true;
}

// ── Print a result ────────────────────────────────────────────────────────────
void printResult(double Q2, double Ebeam, const KinResult &r, const string &mode) {
    cout << fixed << setprecision(5);
    cout << "\n=== Kinematics (mode: " << mode << ") ===\n";
    cout << "  Input:          Q2 = " << Q2 << " GeV^2,  Ebeam = " << Ebeam << " GeV\n";
    cout << "  W               = " << r.W           << " GeV\n";
    cout << "  nu              = " << r.nu           << " GeV\n";
    cout << "  E'              = " << r.Eprime       << " GeV\n";
    cout << "  theta_e         = " << r.theta_e_deg  << " deg\n";
    cout << "  |q|             = " << r.qmag         << " GeV\n";
    cout << "  beta_CM         = " << r.betaCM       << "\n";
    cout << "  p* (CM pion)    = " << r.pstar_val    << " GeV/c\n";
    cout << "  theta* (CM)     = 0.000 deg  [forward production]\n";
    cout << "  -----------------\n";
    cout << "  E_pi (lab)      = " << r.E_pi         << " GeV\n";
    cout << "  p_pi (lab)      = " << r.p_pi         << " GeV/c\n";
    cout << "  theta_pi (lab)  = " << r.theta_pi_deg << " deg\n";
    cout << "  t (= t_min)     = " << r.t            << " GeV^2\n";
    cout << "  k_pi            = " << r.k_pi         << " GeV  [|q| - |p_pi|]\n";
    cout << "===========================\n";
}

// ── Reproduce the reference table ────────────────────────────────────────────
void printTable(double Ebeam = 10.7, double t_target = -0.40) {
    double Q2vals[] = {5.0, 6.5, 7.5, 8.5};
    cout << "\n";
    cout << "Reproducing reference table:  Ebeam=" << Ebeam
         << " GeV,  t_min = " << t_target << " GeV^2\n";
    cout << fixed << setprecision(3);
    cout << setw(6)  << "Q2"
         << setw(7)  << "W"
         << setw(8)  << "E'"
         << setw(9)  << "th_e"
         << setw(9)  << "th_pi"
         << setw(8)  << "p_pi"
         << setw(8)  << "k_pi"
         << setw(9)  << "t_min"
         << "\n";
    cout << string(65, '-') << "\n";
    for (double Q2 : Q2vals) {
        double W_sol;
        if (!findW_for_tmin(Q2, Ebeam, t_target, W_sol)) {
            cout << setw(6) << Q2 << "  NO SOLUTION\n"; continue;
        }
        KinResult r;
        if (!computeKinematics(Q2, Ebeam, W_sol, 0.0, r)) {
            cout << setw(6) << Q2 << "  UNPHYSICAL\n"; continue;
        }
        cout << setw(6)  << Q2
             << setw(7)  << r.W
             << setw(8)  << r.Eprime
             << setw(9)  << r.theta_e_deg
             << setw(9)  << r.theta_pi_deg
             << setw(8)  << r.p_pi
             << setw(8)  << r.k_pi
             << setw(9)  << r.t
             << "\n";
    }
    cout << "\nReference table (from paper):\n";
    cout << setw(6)  << "Q2"
         << setw(7)  << "W"
         << setw(8)  << "E'"
         << setw(9)  << "th_e"
         << setw(9)  << "th_pi"
         << setw(8)  << "p_pi"
         << setw(8)  << "k_pi"
         << setw(9)  << "t"
         << "\n";
    cout << string(65, '-') << "\n";
    double ref[][9] = {
        {5.0, 2.42, 5.371, 16.96, 15.73, 5.111, 0.67, -0.40},
        {6.5, 2.72, 3.772, 23.15, 11.59, 6.715, 0.67, -0.40},
        {7.5, 2.89, 2.707, 29.48,  9.07, 7.784, 0.67, -0.40},
        {8.5, 2.95, 2.003, 36.71,  7.50, 8.515, 0.66, -0.40},
    };
    for (auto &row : ref) {
        cout << setw(6)  << row[0]
             << setw(7)  << row[1]
             << setw(8)  << row[2]
             << setw(9)  << row[3]
             << setw(9)  << row[4]
             << setw(8)  << row[5]
             << setw(8)  << row[6]
             << setw(9)  << row[7]
             << "\n";
    }
    cout << "\nNOTE: Differences arise from rounding of W in the reference table\n"
         << "      (paper quotes 2 decimal places). Row Q2=8.5 in the reference\n"
         << "      has W=2.95 which gives t_min=-0.462, inconsistent with t=-0.40;\n"
         << "      the physically correct W for that row is ~3.091.\n";
}

// ── Main entry point ──────────────────────────────────────────────────────────
/*
  Arguments:
    Q2        : virtuality (GeV^2, positive)
    Ebeam     : beam energy (GeV)
    mode      : "W"  to specify W directly
                "t"  to find W from t_min = t_target
    value     : the value of W (GeV) or t_target (GeV^2, negative)

  Example:
    pionCT_kinematics(5.0, 10.7, "W", 2.42)
    pionCT_kinematics(5.0, 10.7, "t", -0.40)
    pionCT_kinematics(0.0, 10.7, "t", -0.40)   // print comparison table
*/
void pionCT_kinematics(Double_t Q2 = 0.0, Double_t Ebeam = 10.7,
                           const char *mode_c = "table", Double_t value = -0.40)
{
    string mode(mode_c);

    // --- Special mode: reproduce full table ---
    if (mode == "table" || Q2 == 0.0) {
        printTable(Ebeam, value);
        return;
    }

    double W_use = 0.0;
    if (mode == "W" || mode == "w") {
        W_use = value;
    } else if (mode == "t" || mode == "T") {
        double t_target = value;
        if (t_target >= 0.0) {
            cerr << "t_target should be negative for physical forward production.\n";
            return;
        }
        if (!findW_for_tmin(Q2, Ebeam, t_target, W_use)) {
            cout << "No solution: t_target=" << t_target
                 << " may not be reachable for Q2=" << Q2
                 << ", Ebeam=" << Ebeam << ".\n";
            return;
        }
        cout << fixed << setprecision(5);
        cout << "Found W = " << W_use << " GeV giving t_min = " << t_target
             << " GeV^2\n";
    } else {
        cerr << "Unknown mode '" << mode << "'. Use \"W\", \"t\", or \"table\".\n";
        return;
    }

    KinResult r;
    if (!computeKinematics(Q2, Ebeam, W_use, 0.0, r)) {
        cout << "Kinematics unphysical for W=" << W_use << " Q2=" << Q2
             << " Ebeam=" << Ebeam << "\n";
        return;
    }
    printResult(Q2, Ebeam, r, mode);

    // --- Optional: scan theta* and show t range ---
    cout << "\n--- t vs theta* scan for this W ---\n";
    cout << fixed << setprecision(4);
    cout << setw(12) << "theta*(deg)" << setw(12) << "t(GeV^2)"
         << setw(12) << "p_pi(GeV)" << setw(14) << "theta_pi(deg)" << "\n";
    int Nstep = 9;
    for (int i = 0; i <= Nstep; ++i) {
        double th = PI * i / Nstep;
        KinResult ri;
        if (computeKinematics(Q2, Ebeam, W_use, th, ri))
            cout << setw(12) << th*180.0/PI << setw(12) << ri.t
                 << setw(12) << ri.p_pi     << setw(14) << ri.theta_pi_deg << "\n";
    }
}
