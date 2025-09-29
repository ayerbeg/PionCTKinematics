/*
  pionCT_kinematics.C

  Usage examples (from a shell):
    root -l -q "pionCT_kinematics.C(5.0,11.0,-0.40)"
    root -l
    .L pionCT_kinematics.C+
    pionCT_kinematics(5.0,11.0,-0.40);

  Inputs:
   - Q2       : virtuality (GeV^2, positive number)
   - Ebeam    : beam energy (GeV)
   - t_target : chosen Mandelstam t (GeV^2), typically negative (e.g. -0.40)

  Output (printed):
   - W, E', theta_e (deg), nu, |q|, p* (CM pion momentum), p_pi_lab, theta_pi_lab (deg), t_found
*/

#include "TLorentzVector.h"
#include "TMath.h"
#include <iostream>
#include <vector>
#include <limits>
using namespace std;

// Physical constants (GeV units)
const double M_p  = 0.9382720813;   // proton mass
const double M_n  = 0.9395654133;   // neutron mass (not used for constraints here)
const double M_pi = 0.139570611;    // charged pion mass
const double PI   = TMath::Pi();

// Two-body CM momentum for p + pi final state
double pstar_from_W(double W) {
  double s = W*W;
  double term = (s - (M_p + M_pi)*(M_p + M_pi)) * (s - (M_p - M_pi)*(M_p - M_pi));
  if (term <= 0) return 0.0;
  return 0.5/W * sqrt(term);
}

// For given Ebeam, Q2, Eprime and thetaStar (CM), compute pionLab 4-vector and t
// Returns t (GeV^2). Outputs pionLab and thetaPiLab (rad)
double compute_t_for(double Energy_t, double Q2, double Eprime, double thetaStar,
                     TLorentzVector &pionLabOut, double &thetaPiLabOut) {
  return 0.0; // placeholder - will be replaced below
}

// We'll implement compute_t_for below in full (to avoid forward-declare problems)
double compute_t_and_pionLab(double Ebeam, double Q2, double Eprime, double thetaStar,
                             TLorentzVector &pionLabOut, double &thetaPiLabOut) {
  // Basic kinematics
  double nu = Ebeam - Eprime;
  if (nu <= 0) return 1e6;
  double W2 = M_p*M_p + 2.0*M_p*nu - Q2;
  if (W2 <= (M_p + M_pi)*(M_p + M_pi)) return 1e6; // below threshold -> invalid
  double W = sqrt(W2);

  // CM pion momentum
  double pstar = pstar_from_W(W);
  if (pstar <= 0) return 1e6;
  double Epi_star = sqrt(pstar*pstar + M_pi*M_pi);

  // virtual photon 3-momentum magnitude
  double qmag = sqrt(nu*nu + Q2);

  // hadronic CM velocity along q direction (lab frame)
  double betaCM = qmag / (M_p + nu);
  if (betaCM >= 1.0) return 1e6; // unphysical

  // Pion 4-vector in CM (choose phi=0)
  double px_star = pstar * sin(thetaStar);
  double pz_star = pstar * cos(thetaStar);
  TLorentzVector pionStar(px_star, 0.0, pz_star, Epi_star);

  // Boost from CM to lab along +z by betaCM
  TLorentzVector pionLab = pionStar;
  pionLab.Boost(0.0, 0.0, betaCM);


  // virtual photon 4-vector in lab (aligned with +z)
  TLorentzVector qlab(0.0, 0.0, qmag, nu);

  // At this point we've treated the lab z-axis as aligned with the virtual photon (q).
  // In the real lab the beam (initial electron) defines +z and the virtual photon is at
  // the electron scattering angle theta_e relative to the beam. We must rotate the
  // q-aligned "lab" vectors into the true lab frame.

  // compute electron scattering angle (theta_e) from Q2 and energies
  double sin2half = Q2 / (4.0 * Ebeam * Eprime); // use Ebeam parameter here
  double thetaE = 0.0;
  if (sin2half > 0.0 && sin2half <= 1.0) thetaE = 2.0 * asin(sqrt(sin2half));

  // q vector components in the true lab where beam is along +z
  // take initial electron k = (0,0,Ebeam,Ebeam) and scattered k' = (E' sin th, 0, E' cos th, E')
  double kx_prime = Eprime * sin(thetaE);
  double kz_prime = Eprime * cos(thetaE);
  double qx_lab = - kx_prime;                      // qx = kx - kx'
  double qz_lab = Ebeam - kz_prime;             // qz = kz - kz'

  // rotation angle about y that maps +z onto q-hat: alpha = atan2(qx, qz)
  double alpha = atan2(qx_lab, qz_lab);

  // rotate both q and pion from the q-aligned frame into the true lab frame
  qlab.RotateY(alpha);
  pionLab.RotateY(alpha);

  // t = (q - p_pi)^2 computed now in the true lab coordinates
  TLorentzVector diff = qlab - pionLab;
  double t = diff.M2(); // (q - p)^2

// safer / clearer ways to get lab pion kinematics
double ppi_lab = pionLab.P();          // magnitude of lab momentum
double pperp  = pionLab.Pt();          // transverse momentum (Pt) in GeV/c
double thetaPiLab = pionLab.Theta();   // polar angle in radians (0..pi)



  pionLabOut = pionLab;
  thetaPiLabOut = thetaPiLab;
 
  return t;
}

// Try to find theta* (in [0,pi]) that gives t_target for fixed Eprime.
// Use coarse scan to bracket and then bisection refinement.
bool find_thetaStar_for_t_given_Eprime(double Ebeam, double Q2, double Eprime, double t_target,
                                       double &thetaStar_out, TLorentzVector &pionLab_out, double &thetaPiLab_out) {
  const int Nscan = 360;                // coarse grid
  double thmin = 1e-6, thmax = PI-1e-6;
  double prev_th = thmin;
  TLorentzVector prevPion; double prevThetaPi;
  double prev_t = compute_t_and_pionLab(Ebeam,Q2,Eprime,prev_th,prevPion,prevThetaPi);
  if (!isfinite(prev_t)) prev_t = 1e9;
  for (int i = 1; i <= Nscan; ++i) {
    double th = thmin + (thmax - thmin) * i / double(Nscan);
    TLorentzVector pion; double thetaPi;
    double tval = compute_t_and_pionLab(Ebeam,Q2,Eprime,th,pion,thetaPi);
    if (!isfinite(tval)) { prev_th = th; prev_t = tval; continue; }
    // Check if sign change in f(th) = tval - t_target
    double fprev = prev_t - t_target;
    double fcur  = tval - t_target;
    if (fprev * fcur <= 0) {
      // bracket found between prev_th and th -> bisection
      double a = prev_th, b = th;
      double fa = fprev, fb = fcur;
      for (int iter=0; iter<60; ++iter) {
        double c = 0.5*(a+b);
        TLorentzVector pionc; double thpi;
        double fc_t = compute_t_and_pionLab(Ebeam,Q2,Eprime,c,pionc,thpi) - t_target;
        if (!isfinite(fc_t)) break;
        if (fabs(fc_t) < 1e-6) {
          thetaStar_out = c; pionLab_out = pionc; thetaPiLab_out = thpi; return true;
        }
        if (fa * fc_t <= 0) { b = c; fb = fc_t; } else { a = c; fa = fc_t; }
      }
      // fallback: take midpoint
      double cmid = 0.5*(a+b);
      TLorentzVector pionc; double thpi;
      double tmid = compute_t_and_pionLab(Ebeam,Q2,Eprime,cmid,pionc,thpi);
      thetaStar_out = cmid; pionLab_out = pionc; thetaPiLab_out = thpi;
      return true;
    }
    prev_th = th;
    prev_t = tval;
  }
  // no bracket found => return false
  return false;
}

// Main routine: scan E' and find solution (forward-most) matching t_target
void pionCT_kinematics(Double_t Q2, Double_t Ebeam, Double_t t_target) {
  cout << fixed;
  cout.precision(5);

  // E' scanning range and step (adjust if you want more precision)
  double Emin = 0.05;
  double Emax = Ebeam - 0.05;
  double dE = 0.0025; // 2.5 MeV steps; smaller -> slower but more accurate

  bool foundAny = false;
  struct Sol { double Eprime, thetaStar, t, p_pi_lab, theta_pi_lab, W, pstar; TLorentzVector pionLab;};
  vector<Sol> solutions;

  for (double Eprime = Emin; Eprime <= Emax; Eprime += dE) {
    // quick check: Q2 constraint for electron scattering angle (sin^2(theta/2) <=1)
    double denom = 4.0 * Ebeam * Eprime;
    if (denom <= 0) continue;
    double sin2half = Q2 / denom;
    if (sin2half <= 0.0 || sin2half > 1.0) continue; // impossible Eprime for given Q2

    // compute W^2 and threshold condition
    double nu = Ebeam - Eprime;
    double W2 = M_p*M_p + 2.0*M_p*nu - Q2;
    if (W2 <= (M_p + M_pi)*(M_p + M_pi)) continue; // below pi production threshold

    // try to find thetaStar that yields t_target for this Eprime
    double thetaStar_sol; TLorentzVector pionLab_sol; double thetaPiLab_sol;
    bool ok = find_thetaStar_for_t_given_Eprime(Ebeam,Q2,Eprime,t_target,thetaStar_sol,pionLab_sol,thetaPiLab_sol);
    if (ok) {
      double W = sqrt(W2);
      double pstar = pstar_from_W(W);
      TLorentzVector diff; // to recompute exact t
      double t_found = compute_t_and_pionLab(Ebeam,Q2,Eprime,thetaStar_sol,pionLab_sol,thetaPiLab_sol);
      Sol s; s.Eprime = Eprime; s.thetaStar = thetaStar_sol; s.t = t_found;
      s.p_pi_lab = pionLab_sol.P(); s.theta_pi_lab = thetaPiLab_sol;
      s.W = W; s.pstar = pstar; s.pionLab = pionLab_sol;
      solutions.push_back(s);
      foundAny = true;
      // don't break: we want to collect all and choose forward-most later
    }
  }

  if (!foundAny) {
    cout << "No solution found for Q2="<<Q2<<", E="<<Ebeam<<", t="<<t_target<<endl;
    cout << "Try increasing E' scan range, decreasing dE, or relax t tolerance." << endl;
    return;
  }

  // Choose the solution with the smallest thetaStar (forward in CM) -> minimizes FSI
  double bestTheta = 1e9; int bestIdx = -1;
  for (size_t i=0;i<solutions.size();++i) {
    if (solutions[i].thetaStar < bestTheta) { bestTheta = solutions[i].thetaStar; bestIdx = i; }
  }
  Sol best = solutions[bestIdx];

  // compute derived quantities for the chosen solution
  double Eprime = best.Eprime;
  double nu = Ebeam - Eprime;
  double W = best.W;
  double qmag = sqrt(nu*nu + Q2);
  double sin2half = Q2 / (4.0 * Ebeam * Eprime);
  double thetaE = 2.0 * asin(sqrt(sin2half));
  double pstar = best.pstar;





double kpi_CM      = pstar;                         // CM pion momentum (what the macro called k_pi)
double ppi_lab     = best.p_pi_lab;                 // lab momentum magnitude
double ppi_lab_pt  = best.pionLab.Perp();           // lab transverse momentum (p_T)
double theta_pi_lab_deg = best.theta_pi_lab * 180.0 / PI;









  cout << "=== Solution (chosen forward-most) ===\n";
  cout << "Input:  Q2 = " << Q2 << " GeV^2, Ebeam = " << Ebeam << " GeV, t_target = " << t_target << " GeV^2\n";
  cout << "Computed W       = " << W << " GeV\n";
  cout << "Scattered E'     = " << Eprime << " GeV\n";
  cout << "Electron angle   = " << thetaE*180.0/PI << " deg\n";
  cout << "nu (q0)          = " << nu << " GeV, |q| = " << qmag << " GeV\n";

cout << "p*_pi (CM) k_pi      = " << kpi_CM << " GeV/c\n";
cout << "p_pi (lab)           = " << ppi_lab << " GeV/c\n";
cout << "p_pi_transverse (Pt) = " << ppi_lab_pt << " GeV/c\n";
cout << "theta_pi (lab)       = " << theta_pi_lab_deg << " deg\n";

  cout << "t (found)        = " << best.t << " GeV^2\n";
  cout << "theta*_CM (deg)  = " << best.thetaStar*180.0/PI << " deg\n";
  cout << "-------------------------------------\n";
  cout << "Note: If you need more precision, reduce dE in the macro and/or increase Nscan in the theta* scanning function.\n";
}
