class SNInj {
public:
    SNInj();
    ~SNInj();
    
    std::vector<Real> GetSNTimes(Real M, Real t_start, Real t_end, Real unit_time);

private:
    Real a1, a2, a3, a4;
    Real t1, t2, t3;
    Real s1, s2, s3;
    
    Real n_sn(Real t);
    Real get_sn_t(Real n, Real n1, Real n2, Real M);
};

// Constructor
SNInj::SNInj() {
    // From FIRE-3 Paper
    a1 = 0.39;
    a2 = 0.51;
    a3 = 0.18;
    a4 = 0.0083;

    t1 = 3.7; // Myr
    t2 = 7.0;
    t3 = 44.0;

    s1 = std::log(a2/a1)/log(t2/t1);
    s2 = std::log(a3/a2)/std::log(t3/t2);
    s3 = -1.1;
    
    return;
}

// Destructor
SNInj::~SNInj() {
    return;
}

// Number of SN in a given time
Real SNInj::n_sn(Real t) { 
    // Core-Collapse
    if (t < t1) {
        return 0.0;
    } else if (t <= t2) {
        return ((a1*t1)/(s1+1))*(std::pow(t/t1,s1+1)-1)/1000.0;
    } else if (t <= t3) {
        return n_sn(t2) + ((a2*t2)/(s2+1))*(std::pow(t/t2,s2+1)-1)/1000.0;
    }
    
    // Ia
    return n_sn(t3) + ((a4*t3)/(s3+1))*(std::pow(t/t3,s3+1)-1)/1000;
}

// Get the time of the nth SN              
Real SNInj::get_sn_t(Real n, Real n1, Real n2, Real M) {
    if (n < n1) {
        return t1*std::pow(((1000*n*(s1+1))/(M*a1*t1))+1,1/(s1+1));
    } else if (n < n2) {
        return t2*std::pow(((1000*(n-n1)*(s2+1))/(M*a2*t2))+1,1/(s2+1));
    } else {
        return t3*std::pow(((1000*(n-n2)*(s3+1))/(M*a4*t3))+1,1/(s3+1));
    }
    
    return 1e20; // Should never reach here
}

// Get list of SNe. t_start and t_end are in code units. Returns in code units.
std::vector<Real> SNInj::GetSNTimes(Real M, Real t_start, Real t_end, Real unit_time) {
    const Real s_to_Myr = 3.171e-14;

    Real n1 = n_sn(t2)*M;
    Real n2 = n_sn(t3)*M;
    Real n_end = n_sn((t_end-t_start)*unit_time*s_to_Myr)*M;

    std::vector<Real> times;
    for (int i=1; i < n_end; i++) {
        times.push_back(t_start+(get_sn_t(i, n1, n2, M)/(unit_time*s_to_Myr)));
    }
    
    times.push_back(1e20); // Dummy SN at t=inf
    
    return times;
}

// Code below was used for testing - leaving it here for now 

// std::vector<Real> LinearSpacedArray(Real a, Real b, std::size_t N) {
//     Real h = (b - a) / static_cast<double>(N-1);
//     std::vector<Real> xs(N);
//     std::vector<Real>::iterator x;
//     Real val;
//     for (x = xs.begin(), val = a; x != xs.end(); ++x, val += h) {
//         *x = val;
//     }
//     return xs;
// }
    
// // Get a list of supernova times 
// std::vector<Real> get_sn_times(const Real m_sol, const Real t_end) {
//     // Core collaspe
//     std::vector<Real>times = LinearSpacedArray(0.1, 1.1, 1000);
//     times.push_back(1e20);
//     return times;
// }


// Real get_sn_ep(const Real r, const Real nh) {
//     Real E_SN                    = 1.0e51 / energy_scale;
//     Real P_SN                    = 3.479830e42 / mass_scale / vel_scale;
//     Real ejecta_mass             = pin->GetReal("problem", "ejecta_mass") / mass_scale;

//     // Martizzi+15 fitting params
//     Real average_n_H = average_density*rho_scale/(muH*mp);
//     Real alpha       = -7.8 *pow(average_n_H/100.,0.03);
//     Real r_cool      = 3.0  *pow(average_n_H/100.,-0.42) * pc / length_scale;
//     Real r_rise      = 5.5  *pow(average_n_H/100.,-0.40) * pc / length_scale;
//     Real r_0         = 0.97 *pow(average_n_H/100.,-0.33) * pc / length_scale;
//     Real r_break     = 4.0  *pow(average_n_H/100.,-0.43) * pc / length_scale;

//     // Martizzi+15 radial momentum
//     Real P_SN_rad;

//     if ( r < r_break ) {
//         P_SN_rad= P_SN * pow(r/r_0,1.5);
//     } else {
//         P_SN_rad= P_SN * pow(r_break/r_0,1.5);
//     }

//     // Martizzi+15 thermal energy
//     if (r < r_cool ) {
//         E_SN_th = E_SN;
//     } else if (r < r_rise) {
//         E_SN_th = E_SN * pow(r/r_cool,alpha);
//     } else {
//         E_SN_th = E_SN * pow(r_rise/r_cool,alpha);
//     }
    
//     E_SN_th /= (N_cells*vol_cell);
//     //rho_ej = ejecta_mass/vol_cell/N_cells;
//     v_ej   = P_SN_rad/(vol_cell*N_cells*(average_density + rho_ej)) * sqrt((average_density + rho_ej)/average_density);

//     return (E_SN_th,v_ej)
// }
