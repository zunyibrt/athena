std::vector<Real> LinearSpacedArray(Real a, Real b, std::size_t N) {
    Real h = (b - a) / static_cast<double>(N-1);
    std::vector<Real> xs(N);
    std::vector<Real>::iterator x;
    Real val;
    for (x = xs.begin(), val = a; x != xs.end(); ++x, val += h) {
        *x = val;
    }
    return xs;
}
    
// Get a list of supernova times 
std::vector<Real> get_sn_times(const Real m_sol, const Real t_end) {
    // Core collaspe
    std::vector<Real>times = LinearSpacedArray(0.1, 1.1, 100);
    times.push_back(1e20);
    return times;
}

// def get_t_list(M,t_end):
//     // From FIRE-3 Paper
//     Real a1 = 0.39;
//     Real a2 = 0.51;
//     Real a3 = 0.18;
//     Real a4 = 0.0083;

//     Real t1 = 3.7; // Myr
//     Real t2 = 7.0;
//     Real t3 = 44.0;

//     Real s1 = std::log(a2/a1)/std::log(t2/t1);
//     Real s2 = std::log(a3/a2)/std::log(t3/t2);
//     Real s3 = -1.1;

//     Real n1 = n_sn(t2)*M;
//     Real n2 = n_sn(t3)*M;
//     Real n_end = n_sn(t_end)*M;
   
//     def get_sn_t(n):
//         if (n < n1):
//             return t1*np.power(((1000*n*(s1+1))/(M*a1*t1))+1,1/(s1+1))
//         elif (n < n2):
//             return t2*np.power(((1000*(n-n1)*(s2+1))/(M*a2*t2))+1,1/(s2+1))
//         else:
//             return t3*np.power(((1000*(n-n2)*(s3+1))/(M*a4*t3))+1,1/(s3+1))
    
//     last_sn = np.floor(n_end)
    
//     return np.array([get_sn_t(n) for n in np.arange(1,last_sn+1)])


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
