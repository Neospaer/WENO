#include <iostream>
#include <fstream>
#include <cmath>

const int M = 16;
const int N = 10 * M;
const int K_WENO = 3;
const int lo = K_WENO;
const int hi = K_WENO + N - 1;

const double XMIN = 0.;
const double XMAX = 1.;

const double H = (XMAX - XMIN) / N;


double ro[N + 2 * K_WENO];
double ru[N + 2 * K_WENO];
double re[N + 2 * K_WENO];

double ro_old[N + 2 * K_WENO];
double ru_old[N + 2 * K_WENO];
double re_old[N + 2 * K_WENO];

double int_ro[N + 2 * K_WENO];
double int_ru[N + 2 * K_WENO];
double int_re[N + 2 * K_WENO];

double r_[2 * K_WENO];
double p_[2 * K_WENO];
double u_[2 * K_WENO];


const double TMAX = 0.2;

const double TAU = 1.e-3/M;

const double TH = TAU / H;

const double GAM = 1.4;



void WENO(double* UU, double& Um, double& Up);
void calc_flux_hllc(
    double rl, double ul, double pl, double gaml,
    double rr, double ur, double pr, double gamr,
    double& qr, double& qu, double& qe);


int main() {
    /*** Init ***/
    std::cout << H << std::endl;
    for (int i = lo; i <= hi; i++) {
        double x = XMIN + (i - lo + 0.5) * H;
        if (x < 0.5) {
            ro[i] = 1.;
            ru[i] = 0.;
            re[i] = 1. / (GAM - 1.);
        }
        else {
            ro[i] = 0.125;
            ru[i] = 0.;
            re[i] = 0.1 / (GAM - 1.);
        }
    }


    double t = 0.;
    int step = 0;
    while (t < TMAX) {
        step++;
        t += TAU;

        for (int i = lo; i <= hi; i++) {
            ro_old[i] = ro[i];
            ru_old[i] = ru[i];
            re_old[i] = re[i];

            int_ro[i] = 0.;
            int_ru[i] = 0.;
            int_re[i] = 0.;
        }


        /*** Boundary conditions ***/
        for (int i = 0; i < K_WENO; i++) {
            ro[lo - i - 1] = ro[lo + i];
            ru[lo - i - 1] = ru[lo + i];
            re[lo - i - 1] = re[lo + i];

            ro[hi + i + 1] = ro[hi - i];
            ru[hi + i + 1] = ru[hi - i];
            re[hi + i + 1] = re[hi - i];
        }

        /*** Fluxes ***/
        for (int i = lo - 1; i <= hi; i++) {
            double rl, pl, ul;
            double rr, pr, ur;
            double fr, fu, fe;
            for (int k = -K_WENO + 1; k <= K_WENO; k++)
            {
                r_[k + K_WENO - 1] = ro[i + k];
                u_[k + K_WENO - 1] = ru[i + k] / ro[i + k];
                p_[k + K_WENO - 1] = (re[i + k] - 0.5 * ru[i + k] * ru[i + k] / ro[i + k]) * (GAM - 1.);
            }
            WENO(r_, rl, rr);
            WENO(p_, pl, pr);
            WENO(u_, ul, ur);

            calc_flux_hllc(rl, ul, pl, GAM,
                rr, ur, pr, GAM,
                fr, fu, fe);

            int_ro[i] -= fr;
            int_ru[i] -= fu;
            int_re[i] -= fe;

            int_ro[i + 1] += fr;
            int_ru[i + 1] += fu;
            int_re[i + 1] += fe;
        }


        /*** Update fields ***/
        for (int i = lo; i <= hi; i++) {
            ro[i] += TH * int_ro[i];
            ru[i] += TH * int_ru[i];
            re[i] += TH * int_re[i];
        }

        if (step % 100 == 0) {
            std::cout << "STEP = " << step << std::endl;
            std::ofstream fout("WENO160(1).txt");

            fout << "x,r,u,p,e" << std::endl;
            for (int i = lo; i <= hi; i++) {
                double x = XMIN + (i - lo + 0.5) * H;
                double r = ro[i];
                double u = ru[i] / ro[i];
                double e = (re[i] / ro[i] - 0.5 * u * u);
                double p = r * e * (GAM - 1.);
                fout << x << ", " << r << ", " << u << ", " << p << ", " << e << std::endl;
            }

            fout.close();
        }

    }

    return 0;
}




void WENO(double* UU, double& Um, double& Up)
{
    double BETA[3];
    double ALPHA[3];
    double eps = 1.0e-6;
    if ((UU[2] - UU[1]) * (UU[3] - UU[2]) < 0.0) Um = UU[2];
    else
    {
        BETA[0] = (13. / 12.) * (UU[2] - 2 * UU[3] + UU[4]) * (UU[2] - 2 * UU[3] + UU[4]) + 0.25 * (3 * UU[2] - 4 * UU[3] + UU[4]) * (3 * UU[2] - 4 * UU[3] + UU[4]);
        BETA[1] = (13. / 12.) * (UU[1] - 2 * UU[2] + UU[3]) * (UU[1] - 2 * UU[2] + UU[3]) + 0.25 * (UU[1] - UU[3]) * (UU[1] - UU[3]);
        BETA[2] = (13. / 12.) * (UU[0] - 2 * UU[1] + UU[2]) * (UU[0] - 2 * UU[1] + UU[2]) + 0.25 * (UU[0] - 4 * UU[1] + 3 * UU[2]) * (UU[0] - 4 * UU[1] + 3 * UU[2]);
        ALPHA[0] = 0.3 / ((eps + BETA[0]) * (eps + BETA[0]));
        ALPHA[1] = 0.6 / ((eps + BETA[1]) * (eps + BETA[1]));
        ALPHA[2] = 0.1 / ((eps + BETA[2]) * (eps + BETA[2]));
        Um = (ALPHA[0] * (2 * UU[2] + 5 * UU[3] - UU[4]) + ALPHA[1] * (-UU[1] + 5 * UU[2] + 2 * UU[3]) +
            ALPHA[2] * (2 * UU[0] - 7 * UU[1] + 11 * UU[2])) / ((ALPHA[0] + ALPHA[1] + ALPHA[2]) * 6);
    }
    if ((UU[3] - UU[2]) * (UU[4] - UU[3]) < 0.0) Up = UU[3];
    else
    {
        BETA[0] = (13. / 12.) * (UU[3] - 2 * UU[4] + UU[5]) * (UU[3] - 2 * UU[4] + UU[5]) + 0.25 * (3 * UU[3] - 4 * UU[4] + UU[5]) * (3 * UU[3] - 4 * UU[4] + UU[5]);
        BETA[1] = (13. / 12.) * (UU[2] - 2 * UU[3] + UU[4]) * (UU[2] - 2 * UU[3] + UU[4]) + 0.25 * (UU[2] - UU[4]) * (UU[2] - UU[4]);
        BETA[2] = (13. / 12.) * (UU[1] - 2 * UU[2] + UU[3]) * (UU[1] - 2 * UU[2] + UU[3]) + 0.25 * (UU[1] - 4 * UU[2] + 3 * UU[3]) * (UU[1] - 4 * UU[2] + 3 * UU[3]);
        ALPHA[0] = 0.1 / ((eps + BETA[0]) * (eps + BETA[0]));
        ALPHA[1] = 0.6 / ((eps + BETA[1]) * (eps + BETA[1]));
        ALPHA[2] = 0.3 / ((eps + BETA[2]) * (eps + BETA[2]));
        Up = (ALPHA[0] * (11 * UU[3] - 7 * UU[4] + 2 * UU[5]) + ALPHA[1] * (2 * UU[2] + 5 * UU[3] - UU[4]) +
            ALPHA[2] * (-UU[1] + 5 * UU[2] + 2 * UU[3])) / ((ALPHA[0] + ALPHA[1] + ALPHA[2]) * 6);
    }
}




#define F_HLLC_U(UK, FK, SK, SS, PK, RK, VK) (((SS)*((SK)*(UK)-(FK)) + (SK)*( (PK)+(RK)*((SK)-(VK))*((SS)-(VK)) )) / ((SK)-(SS)))
#define F_HLLC_R(UK, FK, SK, SS, PK, RK, VK) (((SS)*((SK)*(UK)-(FK))) / ((SK)-(SS)))
#define F_HLLC_E(UK, FK, SK, SS, PK, RK, VK) (((SS)*((SK)*(UK)-(FK)) + (SK)*( (PK)+(RK)*((SK)-(VK))*((SS)-(VK)) )*(SS)) / ((SK)-(SS)))


void calc_flux_hllc(
    double rl, double ul, double pl, double gaml,
    double rr, double ur, double pr, double gamr,
    double& qr, double& qu, double& qe)
{
    double          sl, sr, p_star, s_star, p_pvrs, ql_, qr_, tmp;

    double czl = sqrt(pl * gaml / rl);
    double czr = sqrt(pr * gamr / rr);

    double el = pl / (rl * (gaml - 1.));
    double er = pr / (rr * (gamr - 1.));

    double e_tot_l = el + 0.5 * (ul * ul);
    double e_tot_r = er + 0.5 * (ur * ur);

    p_pvrs = 0.5 * (pl + pr) - 0.5 * (ur - ul) * 0.25 * (rl + rr) * (czl + czr);
    p_star = (p_pvrs > 0.) ? p_pvrs : 0.;

    ql_ = (p_star <= pl) ? 1 : sqrt(1. + (gaml + 1.) * (p_star / pl - 1.) / (2. * gaml));
    qr_ = (p_star <= pr) ? 1 : sqrt(1. + (gamr + 1.) * (p_star / pr - 1.) / (2. * gamr));

    sl = ul - czl * ql_;
    sr = ur + czr * qr_;

    if (sl > sr) {
        tmp = sl;
        sl = sr;
        sr = tmp;
    }

    s_star = pr - pl;
    s_star += rl * ul * (sl - ul);
    s_star -= rr * ur * (sr - ur);
    s_star /= (rl * (sl - ul) - rr * (sr - ur));

    if (s_star < sl) s_star = sl;
    if (s_star > sr) s_star = sr;


    if (!((sl <= s_star) && (s_star <= sr))) {
        printf("HLLC: inequaluty SL <= S* <= SR is FALSE.\n");
        exit(1);
    }

    if (sl >= 0.) {
        qu = rl * ul * ul + pl;
        qe = (rl * e_tot_l + pl) * ul;
        qr = rl * ul;
    }
    else if (sr <= 0.) {
        qu = rr * ur * ur + pr;
        qe = (rr * e_tot_r + pr) * ur;
        qr = rr * ur;
    }
    else {
        if (s_star >= 0) {
            qu = F_HLLC_U( /*  UK, FK, SK, SS, PK, RK, VK */
                rl * ul,
                rl * ul * ul + pl,
                sl, s_star, pl, rl, ul
            );
            qe = F_HLLC_E( /*  UK, FK, SK, SS, PK, RK, VK */
                rl * e_tot_l,
                (rl * e_tot_l + pl) * ul,
                sl, s_star, pl, rl, ul
            );
            qr = F_HLLC_R( /*  UK, FK, SK, SS, PK, RK, VK */
                rl,
                rl * ul,
                sl, s_star, pl, rl, ul
            );
        }
        else {
            qu = F_HLLC_U( /*  UK, FK, SK, SS, PK, RK, VK */
                rr * ur,
                rr * ur * ur + pr,
                sr, s_star, pr, rr, ur
            );
            qe = F_HLLC_E( /*  UK, FK, SK, SS, PK, RK, VK */
                rr * e_tot_r,
                (rr * e_tot_r + pr) * ur,
                sr, s_star, pr, rr, ur
            );
            qr = F_HLLC_R( /*  UK, FK, SK, SS, PK, RK, VK */
                rr,
                rr * ur,
                sr, s_star, pr, rr, ur
            );
        }
    }
}