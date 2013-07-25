#include "GfcEquilSolver.h"

using namespace std;
using namespace Mutation::Numerics;

namespace Mutation {
    namespace Thermodynamics {

//#define OUTPUT_CONVERGENCE_STATISTICS
#ifdef OUTPUT_CONVERGENCE_STATISTICS
static int newton_count = 0;
static int step_count   = 0;
#endif

// Provide values of the tuneable variables
const double GfcEquilSolver::sm_max_min_frac = 0.01;
const double GfcEquilSolver::sm_logy_lim     = 120.0;
const double GfcEquilSolver::sm_res_tol      = 1.0E-10;
const double GfcEquilSolver::sm_ds_inc       = 3.0;
const double GfcEquilSolver::sm_ds_dec       = 0.25;
const double GfcEquilSolver::sm_dec_min      = 0.75;
const double GfcEquilSolver::sm_err_huge     = 1.0E6;
const double GfcEquilSolver::sm_eps_e        = 1.0E-99;
const double GfcEquilSolver::sm_eps_s        = 1.0E-99;

// Constructor
GfcEquilSolver::GfcEquilSolver(const Thermodynamics& thermo)
    : m_thermo(thermo), m_ns(thermo.nSpecies()), m_ne(thermo.nElements()),
      mp_ls_Br(NULL)
{ 
    int i, k;
    
    // No need to do anything for a single species
    if (m_ns == 1)
        return;
    
    // Form the basic constraint matrix (we only care about elemental
    // constraints right now)
    m_B = m_thermo.elementMatrix();
    
    // Now form the SVD of the basic constraint matrix which can tell us the 
    // the rank of B and which species are determined/undetermined by the 
    // constraints alone
    SVD<double> svd_B(m_B);
    m_nb = svd_B.rank();    // rank of B (ie number of lin-indep columns of B)
    int nc = m_ne;          // total number of columns in B (ie constraints)
    
    // Find determined species (B' is orthogonal to species k direction)
    m_species_order.resize(m_ns, 0);
    m_nsd = 0;
    
    for (i = 0; i < m_ns; ++i) {
        double sum = 0.0;
        for (k = m_nb; k < m_ns; ++k)
            sum += abs(svd_B.U()(i,k));
            
        if (sum < NumConst<double>::eps * svd_B.norm())
            m_species_order[m_nsd++] = i;
    }    
    
    // Order the remaining undetermined species
    m_nsu = m_ns - m_nsd;
    k = 0;
    for (i = 0; i < m_ns; ++i) {
        if (k < m_nsd && m_species_order[k] == i)
            k++;
        else
            m_species_order[m_nsd + i - k] = i;
    }
    
    // Form DEBg (note that we don't actually have the general constraints Bg 
    // but we should stay close to CEQ)
    int ndeg = m_nsd + m_ne;
    RealMatrix DEBg(m_ns, ndeg);
    DEBg(0,m_ns,m_nsd,ndeg) = m_B(0,m_ns,0,nc); // set [E]
    
    for (i = 0; i < m_nsd; ++i) {
        k = m_species_order[i];
        DEBg(k,i) = NumConst<double>::one;               // set [D]
        DEBg(k,k+1,m_nsd,ndeg) = NumConst<double>::zero; // zero row of [E]
    }
    
    // Determine the independent columns of DEBg
    // QR decomposition of DEBg with column pivoting
    QRP<double> qrp(DEBg);
    
    // Figure out the ordering of the columns in AP
    Vector<int> column_order = qrp.columnOrder();
    
    // Determine the independent columns
    Vector<short> independent_columns(ndeg,1);
    for (i = qrp.rank(); i < ndeg; ++i)
        independent_columns(column_order(i)) = 0;
    
    // Identify the determined elements and make them ordered first
    m_element_order.resize(m_ne, 0);
    m_ned = 0;
    
    for (i = 0; i < m_ne; ++i) {
        if (!independent_columns(i+m_nsd))
            m_element_order[m_ned++] = i;
    }
    
    // Order the undetermined elements
    m_neu = m_ne - m_ned;
    k = 0;
    
    for (i = 0; i < m_ne; ++i) {
        if (k < m_ned && m_element_order[k] == i)
            k++;
        else
            m_element_order[m_ned + i - k] = i;
    }
    
    // Assemble the reduced constraint matrix
    m_nrc = m_neu;
    m_B_r = RealMatrix(m_nsu, m_nrc);
    
    for (i = 0; i < m_nsu; ++i) {
        for (k = 0; k < m_neu; ++k) {
            m_B_r(i,k) = 
                m_B(m_species_order[i+m_nsd], m_element_order[k+m_ned]);
        }
    }
    
    // Compute the transpose
    m_B_rt = m_B_r.transpose();
    
    // Store the least-squares object of Br
    mp_ls_Br = new LeastSquares<double>(m_B_r);
    
    // Now compute the modified constraint transformation matrix
    // First reorder U from the SVD of B
    RealMatrix PU(m_ns, m_nb);
    for (i = 0; i < m_ns; ++i) {
        k = m_species_order[i];
        PU(i,i+1,0,m_nb) = svd_B.U()(k,k+1,0,m_nb);
    }
    
    RealMatrix BTPU(m_nb, nc);
    if (m_nsd > 0)
        BTPU(0,m_nsd,0,m_nb) = PU(0,m_nsd,0,m_nb);
    if (m_nsu > 0)
        BTPU(m_nsd,m_nb,0,m_nb) = m_B_rt * PU(m_nsd,m_ns,0,m_nb);

    for (i = 0; i < m_nb; ++i)
        BTPU(0,m_nb,i,i+1) /= svd_B.singularValues()(i);
    
    //m_A = RealMatrix(m_nb, nc);
    m_A = BTPU * svd_B.V().transpose();
    
    // Zero small entries of A
    for (i = 0; i < m_nb; ++i) {
        for (k = 0; k < nc; ++k)
            if (abs(m_A(i,k)) < NumConst<double>::eps * svd_B.norm())
                m_A(i,k) = 0.0;
    }
    
    //cout << "A" << endl;
    //cout << m_A << endl;
    
    // Compute number of atoms in each molecule of undetermined species
    m_au = m_B_r * RealVector(m_nrc, 1.0);
}


pair<int, int> GfcEquilSolver::equilibrate(
    const double &T, const double &P, const double *const p_ev, 
    double *const p_sv)
{

#ifdef OUTPUT_CONVERGENCE_STATISTICS
    newton_count = 0;
    step_count   = 0;
#endif
    pair<int,int> stats = make_pair(0, 0);
    
    // Special case for a single species mixture
    if (m_ns == 1) {
        p_sv[0] = p_ev[0];
        return make_pair(0, 0);
    }
    
    static RealVector cb(m_ne);
    static RealVector cm(m_nb);
    static RealVector zd(m_nsd);
    static RealVector cr(m_nrc);
    static RealVector nmm(m_nsu);
    static RealVector Q(m_nsu);
    
    //for (int i = 0; i < m_ne; ++i)
    //    cout << i << " " << p_ev[i] << endl;
    
    // Initialize constraint dependent quantities if the constraints are
    // different from the last call
    if (cb != asVector(p_ev, m_ne)) {
        // Basic constraints
        cb = asVector(p_ev, m_ne);
        
        // Modified constraint vector
        cm = m_A * asVector(p_ev, m_ne);
        
        // Determined species
        if (m_nsd > 0)
            zd = cm(0,m_nsd);
        
        // Reduced constraint vector
        cr = cm(m_nsd,m_nb);
        
        // Possibly perturb the reduced constraint vector and get the max-min
        // composition of the undetermined species
        perturb(cr, nmm);
        
        // Compute the normalization condition vector Q (Q'Xu = 1)
        Q = m_au * (zd.sum() / cr.sum()) + 1.0;
    }
    
    // Compute the Min-g solution
    RealVector ng(m_nsu);
    RealVector gu(m_nsu);
    ming(T, P, cr, gu, ng);
    
    // Initial solution for undetermined species combines the Min-g and max-min
    // solutions via a linear combination
    RealVector zu0 = ng + sm_max_min_frac * (nmm - ng);
    
    // Initialize gu0 and lam0
    RealVector lam0(m_nrc);
    RealVector gu0(m_nsu);
    lamgInit(zd, zu0, gu, lam0, gu0);
    
    RealVector zu(m_nsu);
    solveFixedT(gu, lam0, gu0, 0.0, Q, cr, zu, stats);
    
    // Finally, reorder the species moles to their original order
    for (int i = 0; i < m_nsd; ++i)
        p_sv[m_species_order[i]] = zd(i);
    for (int i = 0; i < m_nsu; ++i)
        p_sv[m_species_order[i+m_nsd]] = zu(i);

    //cout << "Before returning:" << endl;
    //for (int i = 0; i < m_ns; ++i)
    //    cout << setw(14) << p_sv[i];
    //cout << endl;

    return stats;
} 

void GfcEquilSolver::solveFixedT(
    const RealVector& gu, const RealVector& lam0, const RealVector& gu0, 
    const double& res0, const RealVector& Q, const RealVector& cr,
    RealVector& zu, pair<int,int>& stats) const
{
    // Initialize prior to time steps
    double s = 0.0;
    double send = 1.0;
    RealVector lam(lam0);
    double irr_error = 0.0;
    RealVector dguds(gu);
    dguds -= gu0;
    double error_tol = max( res0 + sm_res_tol, 1.05 * res0 );
    double error = 0.0;
    
    double ds;
    RealVector gus(m_nsu);
    RealVector dlamds(m_nrc);
    RealVector lam1(m_nrc);
    RealVector lam2(m_nrc);
    RealVector y(m_nsu);
    
    // Integrate from s = 0 to s = 1
    while (s < 1.0) {
        ds = send - s;

#ifdef OUTPUT_CONVERGENCE_STATISTICS        
        step_count++;
        cout << "time: s, ds, steps, newts, error = " 
             << s << "  " << ds << "  " << step_count << "  " << newton_count
             << "  " << error << endl;
#endif
        stats.first++;
        
        // Take explicit Euler step (from s to end)
        gus = gu0 + s * dguds;
        
        if (!computeRate(lam, gus, dguds, Q, dlamds)) {
            // Overflow in evaluating y (repeat step with smaller ds)
            send = s + ds * sm_ds_dec;
            continue;
        }
        
        lam1 = lam + ds * dlamds;
        gus  = gu0 + send * dguds;
        
        // Perform Newton iteration to convergence
        if (!newton(cr, gus, lam1, Q, lam2, y, error, stats) || 
            error > error_tol) {
            // Reduce time step and repeat if Newton did not converge
            send = s + ds * sm_ds_dec;
            continue;
        }
        
        // Acceptable convergence: set acceptable tolerance for next step
        error_tol = max(error + sm_res_tol, 1.05 * error);
        
        // Update solution
        lam  = lam2;
        s    = send;
        send = min(send + ds * sm_ds_inc, 1.0);
    } // while (s < 1.0)
    
    // Converged solution obtained s=1
    zu = square(y);    
    RealVector v = m_B_rt * zu;
    zu *= (cr.norm2() / v.norm2());
}

bool GfcEquilSolver::newton(
    const RealVector& cr, const RealVector& gu, const RealVector& lam0, 
    const RealVector& Q, RealVector& lam, RealVector& y, double& error,
    pair<int,int>& stats) const
{
    static RealVector v(m_nrc);
    static RealVector dw(m_nrc);
    static RealVector lam_v(m_nrc);
    static RealVector lam_w(m_nrc);
    static RealVector lam_best(m_nrc);
    static RealVector y_best(m_nsu);
    static RealVector YQ(m_nsu);
    static RealVector P(m_nrc);    
    static RealMatrix H(m_nsu, m_nrc);
    
    const double cr_norm = cr.norm2();
    double err_best, err_last;
    double v_norm;
    double q;
    double v_star, w_star;
    double dalpha;
    
    // Initialize before starting the iterations
    error     = 2.0 * sm_err_huge;
    err_last  = error;
    err_best  = error;  
    lam       = lam0;
    
    int i;
 
    // Iterate until exit due to (a) y failure, (b) error < error_tol, or 
    // (c) no longer converging
    while (true) {
        // Form y
        if (!formy(gu, lam, y)) {
            error = 2.0 * sm_err_huge;
            break; // previous solution may exist
        }
        
        //cout << "y = " << endl << y << endl;

        // Form H
        H = m_B_r;   
        for (i = 0; i < m_nsu; ++i)
            H.row(i) *= y(i);
        //cout << "H = " << endl << H << endl;
        
        // Form v, ||r|| = error, dw, YQ, and q
        v = y * H;
        v_norm = v.norm2();
        //cout << "|v| = " << v_norm << endl;
        dw = v / v_norm - cr / cr_norm;
        //cout << "dw = " << dw << endl;
        
        error = dw.norm2();
        if (error < err_best) {
            err_best = error;
            lam_best = lam;
            y_best = y;
        }
        
        YQ = Q * y;
        //cout << "YQ = " << endl << YQ << endl;      
        q = 1.0 - dot(YQ, y);
        //cout << "q = " << q << endl;

#ifdef OUTPUT_CONVERGENCE_STATISTICS        
        newton_count++;
        cout << "  newt: q, r, newts = " 
             << q << "  " << error << "  " << newton_count << endl;
#endif
        stats.second++;
        
        // Check for convergence (if done, return with current solution)
        if (error + abs(q) < sm_res_tol)
            return true;
        
        // Check to see if we are still making progress
        if (error > sm_dec_min * err_last)
            break;
        
        err_last = error;
        dw *= -v_norm;
        
        // Compute lam_v and lam_w
        LeastSquares<double> ls_H(H);
        ls_H.solve(lam_v, y);
        ls_H.solveATA(lam_w, dw);
        
        // Compute P
        P = YQ * H;
        //cout << "P = " << endl << P << endl;
        
        // Compute P'lam_w and P'lam_v
        w_star = dot(P, lam_w);
        v_star = dot(P, lam_v);
        
        // Now for dalpha
        dalpha = (q - w_star) / (v_star + w_star > 0 ? 
            max(v_star + w_star, 0.5) : min(v_star + w_star, -0.5));
        /*cout << "newton" << endl;
        cout << dalpha << endl;
        cout << q << endl;
        cout << w_star << endl;
        cout << v_star << endl;
        cout << "lam_w = " << endl << lam_w << endl;
        cout << "lam_v = " << endl << lam_v << endl;*/        
        
        // Finally update lam
        lam += (1.0 + dalpha) * lam_w + dalpha * lam_v;
    } // end of Newton iterations
    
    if (err_best > sm_err_huge)
        return false;
    
    // Return with best solution
    if (err_best < error) {
        error = err_best;
        lam = lam_best;
        y = y_best;
    }
    
    return true;
}

void GfcEquilSolver::lamgInit(
    const RealVector& zd, const RealVector& zu0, const RealVector& gu,
    RealVector& lam0, RealVector& gu0) const
{
    // Solve least squares problem min|| -gu + Br*lam0 - ln(xu0)||_2 to get lam0
    RealVector lnxu = log(zu0 / (zd.sum() + zu0.sum()));
    RealVector rhs = gu + lnxu;
    mp_ls_Br->solve(lam0, rhs);
    
    // Finally setup gu0 to satisfy Br*lam0 - ln(xu0) = gu(0)
    gu0 = m_B_r * lam0 - lnxu;
}

void GfcEquilSolver::perturb(RealVector& cr, RealVector& nmm) const
{
    // Possibly perturb the elemental constraints (only constraints for us)
    double ne_max = cr.maxValue();
    double ne_low = sm_eps_e * ne_max;    
    cr = Numerics::max(cr, ne_low);
    
    //cout << ne_max << endl;
    //cout << ne_low << endl; 
    //cout << cr << endl;
    
    // Determine upper bound on moles of undetermined species (ie a species
    // contains all the moles of a given element and thus cannot be any bigger)
    RealVector nk_max(m_nsu, 0.0);
    for (int i = 0; i < m_nsu; ++i) {
        for (int j = 0; j < m_neu; ++j)
            nk_max(i) = std::max(nk_max(i), (m_B_r(i,j) / cr(j)));
        nk_max(i) = 1.0 / nk_max(i);
    }

    // Compute the max-min composition of the undetermined species using the
    // perturbed constraints
    maxmin(cr, nmm);
    
    // Compute a (possibly) perturbed vector of undetermined species moles
    // (just storing in nk_max)
    nk_max = max(nmm, nk_max * sm_eps_s);
    
    // Now back out the (possibly) perturbed constraint vector
    cr = m_B_rt * nk_max;
}

void GfcEquilSolver::maxmin(const RealVector& c, RealVector& nmm) const
{
    // Setup the tableau for input to lp()
    RealMatrix tableau(m_nrc, m_nsu+1);
    
    // B'
    tableau(0,m_nrc,0,m_nsu) = m_B_rt;    
    
    // sum(B)'
    for (size_t i = 0; i < m_nrc; ++i)
        for (size_t j = 0; j < m_nsu; ++j)
            tableau(i,m_nsu) += tableau(i,j);
    
    // Objective function
    RealVector f(m_nsu+1);
    f(m_nsu) = 1.0;
    
    // Temporary solution vector (to make room for zmin)
    RealVector x(m_nsu+1);
    
    // Now solve the LP problem
    double zmin;
    LpResult lp_result = lp(f, MAXIMIZE, tableau, c, 0, 0, x, zmin);
    
    if (lp_result != SOLUTION_FOUND) {
        cout << "Error finding max-min solution!" << endl;
        exit(1);
    }
    
    // Transfer solution
    nmm = x(0,m_nsu) + zmin;
}

void GfcEquilSolver::ming(
    const double& T, const double& P, const RealVector& cr, RealVector& gu,
    RealVector& ng) const
{
    double g[m_ns];
    
    // Compute the standard state normalized molar specific Gibbs energy for
    // each species in the mechanism
    m_thermo.speciesGOverRT(T, P, g);
    
    // Solve the linear programming problem
    LpResult lp_result =
        lp(gu, MINIMIZE, m_B_rt, cr, 0, 0, ng, g[0]);
    
    if (lp_result != SOLUTION_FOUND) {
        cout << "Error finding min-g solution!" << endl;
        exit(1);
    }
}

bool GfcEquilSolver::computeRate(
    const RealVector& lam, const RealVector& gu, const RealVector& dgudt, 
    const RealVector& Q, RealVector& dlamdt) const
{
    // Form the y vector
    RealVector y(m_nsu);
    if (!formy(gu, lam, y))
        return false;
    
    // Form the H matrix (diag(y) * Br), and Y*dg/dt and Q'*Y vectors
    RealMatrix H = m_B_r;    
    for (int i = 0; i < m_nsu; ++i)
        H(i,i+1,0,m_nrc) *= y(i);
    
    // Use Ydgdt before it gets destroyed to compute by the least squares solver    
    RealVector Ydgdt = y * dgudt;
    RealVector QTY = y * Q;
    RealVector QTYH = m_B_rt * (QTY * y);
    double alpha = dot(QTY, Ydgdt);
    
    // Solve the least squares problems for lamdot_g and lamdot_y   
    RealVector lamdot_g(m_nrc);
    RealVector lamdot_y(m_nrc);
    
    LeastSquares<double> ls_H(H);
    ls_H.solve(lamdot_g, Ydgdt);
    ls_H.solve(lamdot_y, y);
    
    // Now finish computing alpha
    alpha = (alpha - dot(QTYH, lamdot_g)) / dot(QTYH, lamdot_y);
    
    // Finally compute the rate
    dlamdt = alpha * lamdot_y + lamdot_g;    
    return true;
}

bool GfcEquilSolver::formy(
    const RealVector& gu, const RealVector& lam, RealVector& y) const
{
    // Form log(y)
    y = 0.5 * (m_B_r * lam - gu);
    
    // Check for overflow
    if (y.maxValue() > sm_logy_lim)
        return false;
    
    // Compute y
    y.exp();
    return true;
}

} // namespace Thermodynamics
} // namespace Mutation


