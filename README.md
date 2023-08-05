# hybridNanofluidOpenLB
This is a code for simulation of hybrid nanofluid written based on OpenLB-1.4.


## Mathematical Relationships
$$ \phi_{np} = \phi_{np1} + \phi_{np2} $$
$$ \rho = \left{1 - \phi_{np}}\right) \rho_{f} + \phi_{np1} \rho_{np1} + \phi_{np2} \rho_{np2} $$
$$ C_p =  {{\left({1 - \phi_{np}}\right) \rho_{f} C_{p, f} + \phi_{np1} \rho_{np1} C_{p,np1} + \phi_{np2} \rho_{np2} C_{p,np2}} \over \rho}  $$
$$ kappa_ = ( ( kappanp1_ + ( lambdanp1_ - 1.0 ) * kappabf_ - ( lambdanp1_ - 1.0 ) * ( kappabf_ - kappanp1_ ) * \phi_{np1} ) / ( kappanp1_ + ( lambdanp1_ - 1.0 ) * kappabf_ + ( kappabf_ - kappanp1_ ) * \phi_{np1} ) ) * ( ( kappanp2_ + ( lambdanp2_ - 1.0 ) * kappabf_ - ( lambdanp2_ - 1.0 ) * ( kappabf_ - kappanp2_ ) * \phi_{np2} ) / ( kappanp2_ + ( lambdanp2_ - 1.0 ) * kappabf_ + ( kappabf_ - kappanp2_ ) * \phi_{np2} ) ) * kappabf_;
$$ mu_ = mubf_ / pow( 1.0 - \phi_{np}, 2.5) $$
$$ beta_ = ((1 - \phi_{np}) * \rho_{f} * betabf_ + \phi_{np1} * \rho_{np1} * betanp1_ + \phi_{np2} * \rho_{np2} * betanp2_)/rho_ $$
$$ nu_ = mu_ / \rho $$
$$ alpha_ = kappa_ / (\rho * C_p) $$
$$ Pr_ = nu_ / alpha_ $$





## Installation
It is working on OpenLB-1.4.
```bash
git clone https://github.com/EhsanGLB/hybridNanofluidOpenLB.git
cd hybridNanofluidOpenLB/hybridNanofluidOpenLB
make
```


## Getting Started
Change the directory of root in definitions.mk file.
```bash
./thermalNaturalConvection2D
```


## References
* [Golab, Ehsan, Behzad Vahedi, Ankur Jain, Robert A. Taylor, and Kambiz Vafai. "Laminar forced convection in a tube with a nano-encapsulated phase change materials: Minimizing exergy losses and maximizing the heat transfer rate." Journal of Energy Storage 65 (2023): 107233.](https://www.sciencedirect.com/science/article/abs/pii/S2352152X23006308)
* [Vahedi, Behzad, Ehsan Golab, Arsalan Nasiri Sadr, and Kambiz Vafai. "Thermal, thermodynamic and exergoeconomic investigation of a parabolic trough collector utilizing nanofluids." Applied Thermal Engineering 206 (2022): 118117.](https://www.sciencedirect.com/science/article/abs/pii/S1359431122000813)
* [Golab, Ehsan, Sahar Goudarzi, Hamed Kazemi-Varnamkhasti, Hossein Amigh, Ferial Ghaemi, Dumitru Baleanu, and Arash Karimipour. "Investigation of the effect of adding nano-encapsulated phase change material to water in natural convection inside a rectangular cavity." Journal of Energy Storage 40 (2021): 102699.](https://www.sciencedirect.com/science/article/abs/pii/S2352152X21004357)
* [Abbasi, Mohammad, Amin Nadimian Esfahani, Ehsan Golab, Omid Golestanian, Nima Ashouri, S. Mohammad Sajadi, Ferial Ghaemi, Dumitru Baleanu, and A. Karimipour. "Effects of Brownian motions and thermophoresis diffusions on the hematocrit and LDL concentration/diameter of pulsatile non-Newtonian blood in abdominal aortic aneurysm." Journal of Non-Newtonian Fluid Mechanics 294 (2021): 104576.](https://www.sciencedirect.com/science/article/abs/pii/S0377025721000859)
* [Jafarzadeh, Sina, Arsalan Nasiri Sadr, Ehsan Kaffash, Sahar Goudarzi, Ehsan Golab, and Arash Karimipour. "The effect of hematocrit and nanoparticles diameter on hemodynamic parameters and drug delivery in abdominal aortic aneurysm with consideration of blood pulsatile flow." Computer Methods and Programs in Biomedicine 195 (2020): 105545.](https://www.sciencedirect.com/science/article/abs/pii/S0169260720307914)
* [Goudarzi, Sahar, Masih Shekaramiz, Alireza Omidvar, Ehsan Golab, Arash Karimipour, and Aliakbar Karimipour. "Nanoparticles migration due to thermophoresis and Brownian motion and its impact on Ag-MgO/Water hybrid nanofluid natural convection." Powder Technology 375 (2020): 493-503.](https://www.sciencedirect.com/science/article/abs/pii/S0032591020307397)
* [Motlagh, Saber Yekani, Ehsan Golab, and Arsalan Nasiri Sadr. "Two-phase modeling of the free convection of nanofluid inside the inclined porous semi-annulus enclosure." International Journal of Mechanical Sciences 164 (2019): 105183.](https://www.sciencedirect.com/science/article/abs/pii/S0020740319315279)
