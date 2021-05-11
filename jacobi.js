class jacobi {
    // Single precision complete elliptic integral of the first kind
    // T. Fukushima, "Fast Computation of Complete Elliptic Integrals
    // and Jacobian Elliptic Functoons", Celest. Mech. Dyn. Astron.,
    // 105, 2009, pp. 305-328.
    //
    // Ported to JS by D.Reznik, May 2021.     
    // Input: mc = complementary parameter 0 <= mc <= 1
    static relk(mc) {
        let P, Q;
        const D1 = 1.0 / 16.0, D2 = 1.0 / 32.0, D3 = 21.0 / 1024.0;
        const D4 = 31.0 / 2048.0, D5 = 6257.0 / 524288.0, D6 = 10293.0 / 1048576.0;
        const PIHALF = 1.57079633, PIINV = 0.318309886;
        let m = 1.0 - mc
        if (mc < 1.05e-8) {
            return 1.38629436 - 0.5 * log(mc);
        }
        if (mc < 0.1) {
            let nome = mc * (D1 + mc * (D2 + mc * (D3 + mc * (D4 + mc * (D5 + mc * D6)))))
            let mx = mc - 0.05;
            let kkc = 1.59100345 + mx * (0.41600074 + mx * (0.24579151 + mx * (0.17948148 + mx * (0.14455606))));
            return -kkc * PIINV * log(nome);
        }
        if (m <= 0.1) {
            let mx = m - 0.05;
            return 1.59100345 + mx * (0.41600074 + mx * (0.24579151 + mx * (0.17948148 + mx * (0.14455606))));
        }
        if (m <= 0.2) {
            let mx = m - 0.15;
            return 1.63525673 + mx * (0.47119063 + mx * (0.30972841 + mx * (0.25220831 + mx * (0.22672562))));
        }
        if (m <= 0.3) {
            let mx = m - 0.25;
            return 1.68575035 + mx * (0.54173185 + mx * (0.40152444 + mx * (0.36964247 + mx * (0.37606072))));
        }
        if (m <= 0.4) {
            let mx = m - 0.35;
            return 1.74435060 + mx * (0.63486428 + mx * (0.53984256 + mx * (0.57189271 + mx * (0.67029514 + mx * (0.83258659)))));
        }
        if (m <= 0.5) {
            let mx = m - 0.45;
            return 1.81388394 + mx * (0.76316325 + mx * (0.76192861 + mx * (0.95107465 + mx * (1.31518067 + mx * (1.92856069)))));
        }
        if (m <= 0.6) {
            let mx = m - 0.55;
            return 1.89892491 + mx * (0.95052179 + mx * (1.15107759 + mx * (1.75023911 + mx * (2.95267681 + mx * (5.28580040)))));
        }
        if (m <= 0.7) {
            let mx = m - 0.65;
            return 2.00759840 + mx * (1.24845723 + mx * (1.92623466 + mx * (3.75128964 + mx * (8.11994455 + mx * (18.6657213 + mx * (44.6039248))))));
        }
        if (m <= 0.8) {
            let mx = m - 0.75;
            return 2.15651565 + mx * (1.79180564 + mx * (3.82675129 + mx * (10.3867247 + mx * (31.4033141 + mx * (100.923704 + mx * (337.326828 + mx * (1158.70793)))))));
        }
        if (m <= 0.85) {
            let mx = m - 0.825;
            return 2.31812262 + mx * (2.61692015 + mx * (7.89793508 + mx * (30.5023972 + mx * (131.486937 + mx * (602.984764 + mx * (2877.02462))))));
        } // all failed
        let mx = m - 0.875;
        return 2.47359617 + mx * (3.72762424 + mx * (15.6073930 + mx * (84.1285084 + mx * (506.981820 + mx * (3252.27706 + mx * (21713.2424 + mx * (149037.045)))))));
    }
    //
    //	Single precision subroutine to compute three Jacobian elliptic functions simultaneously
    //  For limited argument: 0 <= u < K/2
    //
    //    Reference: T. Fukushima, (2012) Numer. Math. DOI 10.1007/s00211-012-0498-0
    //      "Precise and Fast Computation of Jacobian Elliptic Functions by
    //        Conditional Duplication"
    //     Author: T. Fukushima Toshio.Fukushima@nao.ac.jp
    //
    //     Inputs: u = argument, mc = 1-m, 0 < mc <= 1
    //     Output: s = sn(u|m), c=cn(u|m), d=dn(u|m)
    static rscd2(u, mc) {
        //    real*4 u,mc,s,c,d
        //   real*4 m,uA,uT,u0,v,a,b,y,z,my,mc2,m2,xz,w
        //    integer n,j,i
        const m = 1.0 - mc, uA = 1.76269 + mc * 1.16357, uT = 9.207e-4 - m * 4.604e-4;
        let u0 = u;
        for (let n = 0; n <= 20; n++) {
            if (u0 < uT) break;
            u0 = u0 * 0.5
        }
        if (n == 20) console.log("(rscd2) Too large input argument: u=", u);
        let v = u0 * u0, a = 1.0, b = v * 0.5;
        let y, z, my;
        let goto2 = false;
        if (u < uA) {
            for (let j = 1; j <= n; j++) {
                y = b * (a * 2.0 - b); z = a * a; my = m * y; b = (y * 2.0) * (z - my); a = z * z - my * y;
            } else for (let j = 1; j <= n; j++) {
                y = b * (a * 2.0 - b); z = a * a; my = m * y;
                if (z < my * 2.0) { goto2 = true; break; }
                b = (y * 2.0) * (z - my); a = z * z - my * y
            }
        }
        if (!goto2) {
            b = b / a; y = b * (2.0 - b);
            return { c: 1.0 - b, s: Math.sqrt(y); d: Math.sqrt(1.0 - m * y) };
        }
        let c = a - b, mc2 = mc * 2.0, m2 = m * 2.0;
        for (let i = j; i <= n; i++) {
            x = c * c; z = a * a; w = m * x * x - mc * z * z; xz = x * z; c = mc2 * xz + w; a = m2 * xz - w;
        }
        let c = c / a; x = c * c;
        return { c: c, s: Math.sqrt(1.0 - x), d: Math.sqrt(mc + m * x) };
    }
}