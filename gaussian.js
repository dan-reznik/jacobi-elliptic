jacobi_curvature = class {
  constructor(a) {
    this.a = a;
    [this.cA, this.cB] = caustic_axes(a);
    this.d = Math.sqrt(this.cA*this.cA - this.cB*this.cB);
    this.m = this.d/this.cA;
    this.k = jacobi.ellipticK(this.m*this.m);
  }
  // v: (1-v)*p1+v*p2 
  // u: from 0 to 4K
  // ae, be: outer ellipse semi-axes
  // m: caustic eccentricity, stachel
  // K: quarter period
  // i: vtx: 0,1,2
  static gaussian_curvature(v,u,ae,be,m,K,i) { 
    const scd4 = jacobi.sn_cn_dn(u + ((1+i)* 4 * K) / 3, m);
    const scd8 = jacobi.sn_cn_dn(u + ((2+i)* 4 * K) / 3, m);
    const [sn4, cn4, dn4] = [scd4.sn, scd4.cn, scd4.dn];
    const [sn8, cn8, dn8] = [scd8.sn, scd8.cn, scd8.dn];
    const sn82 = sn8 * sn8;
    const sn42 = sn4 * sn4, sn84 = sn82 * sn82, m2 = m * m;
    const ae2 = ae * ae, be2 = be * be, vm12 = (v - 1) * (v - 1), v2 = v * v;
    const sn44=sn42*sn42;
    const sn43=sn4*sn42;
    const delta = -2 * m2 * ae2 * vm12 * (sn82 - 1 / 2) * be2 * sn44 -
          2 * sn8 * m2 * ae2 * be2 * vm12 * (cn8 * cn4 - 1) *
          sn43 + (-2 * sn84 * ae2 * be2 * m2 * v2 +
                  4 * ae2 * (dn8 * v * (-1 + v) * dn4 + ((m2 + 2) * (v2 - v + 1 / 2)) / 2) * be2 *
                  sn82 + 2 * m2 * ae2 * be2 * cn8 * vm12 * cn4 -
                  2 * dn8 * v * ae2 * be2 * (-1 + v) *
                  dn4 + (-1 + (-2 * m2 * v2 + 4 * m2 * v - 2 * m2 - 2 * v2 + 2 * v - 1) * ae2) *
                  be2 + ae2) * sn42 +
          4 * ae2 * (-m2 * v2 * be2 * (cn8 * cn4 - 1) *
                     sn82 / 2 + (dn8 * v * (-1 + v) * dn4 + v2 - v + 1 / 2) * cn8 * be2 * cn4 -
                     dn8 * v * be2 * (-1 + v) * dn4 - 1 / 2 + (-v2 + v - 1 / 2) * be2) * sn8 * sn4 +
          sn84 * ae2 * be2 * m2 *
          v2 + (2 * m2 * v2 * ae2 * be2 * cn4 * cn8 -
                2 * dn8 * v * ae2 * be2 * (-1 + v) *
                dn4 + (-1 + (-1 + (-2 * m2 - 2) * v2 + 2 * v) * ae2) * be2 + ae2) *
          sn82 - 4 * (dn8 * v * ae2 * (-1 + v) * dn4 + 1 / 2 + (v2 - v + 1 / 2) * ae2) *
          be2 * (cn8 * cn4 - 1);
    const f1 = ae * be * (dn4 + dn8) * (sn4 * sn8 + cn8 * cn4 - 1) / delta;
    return -f1 * f1;
  }
  gaussian(v12,t,edge012) {
    const u = 2*this.k*t/Math.PI;
    return this.constructor.gaussian_curvature(v12, u, this.a, 1, this.m, this.k,edge012);
  }
  gaussian3(v12,t) {
    const u = 2*this.k*t/Math.PI;
    return [0,1,2].map(i=>this.constructor.gaussian_curvature(v12, u, this.a, 1, this.m, this.k,i));
  }
  static mean_curvature(v,u,ae,be,m,K,i) { 
    const scd4 = jacobi.sn_cn_dn(u + ((1+i)* 4 * K) / 3, m);
    const scd8 = jacobi.sn_cn_dn(u + ((2+i)* 4 * K) / 3, m);
    const [sn4, cn4, dn4] = [scd4.sn, scd4.cn, scd4.dn];
    const [sn8, cn8, dn8] = [scd8.sn, scd8.cn, scd8.dn];
    const sn42 = sn4*sn4, sn43 = sn4*sn42;
    const sn82 = sn8*sn8, sn84 = sn82*sn82;
    const m2 = m*m, ae2 = ae*ae, be2 = be*be, v2 = v*v; 
    const delta = (-2*((sn82 - 1/2)*sn42 + sn8*(cn4*cn8 - 1)*sn4 - sn82/2 - 
                       cn8*cn4 + 
                       1)*((m2*sn42 + m2*sn82 - 2*dn4*dn8 - 2)*
                           v2 + (-2*m2*sn42 + 2*dn4*dn8 + 2)*v + m2*sn42 - 1)*ae2 - 
                   sn82 - 2*cn8*cn4 - sn42 + 2)*be2 + ae2*(sn8 - sn4)*(sn8-sn4); 
    const numer = (2*ae*
                   be*(dn4 + dn8)*(cn4*cn8 + sn4*sn8 - 
                                   1)*((((ae2 - be2)*sn4 - sn8*ae2)*cn4 + sn4*be2*cn8)*
                                       dn4 - (-sn8*be2*cn4 + (sn4*ae2 + (-ae2 + be2)*sn8)*cn8)*
                                       dn8) - 2*
                   ae*(cn8*m2*sn43 - 
                       sn8*cn4*m2*sn42 + (m2*cn4/2 - (m2*sn82 + 1/2*m2)*cn8)*
                       sn4 + ((m2*sn82 - 1/2*m2)*cn4 + cn8*m2/2)*sn8)*
                   be*((ae2 - be2)*sn42 - 2*ae2*sn8*sn4 + (ae2 - be2)*sn82 - 
                       2*be2*(cn4*cn8 - 1)))*v + 
          2*ae*be*(dn4 + dn8)*(cn4*cn8 + sn4*sn8 - 
                               1)*(-((ae2 - be2)*sn4 - sn8*ae2)*cn4 - sn4*be2*cn8)*dn4 - 
          2*ae*(-cn8*m2*sn43 + 
                sn8*cn4*m2*sn42 + (-m2*cn4/2 - (-1/2 - m2/2)*cn8)*sn4 - 
                sn8*cn4/2)*
          be*((ae2 - be2)*sn42 - 2*ae2*sn8*sn4 + (ae2 - be2)*sn82 - 
              2*be2*(cn4*cn8 - 1));
    return numer/Math.sqrt(delta*delta*delta)
  }
  mean(v12,t,edge012) {
    const u = 2*this.k*t/Math.PI;
    return this.constructor.mean_curvature(v12, u, this.a, 1, this.m, this.k,edge012);
  }
  mean3(v12,t) {
    const u = 2*this.k*t/Math.PI;
    return [0,1,2].map(i=>this.constructor.mean_curvature(v12, u, this.a, 1, this.m, this.k,i));
  }
}