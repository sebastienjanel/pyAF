# Models and formulas of pyAF

Formulas used for each computation, with their references.

Notation: $F$ force, $\delta$ indentation, $E$ Young's modulus, $\nu$ Poisson ratio, $R$ tip radius, $\theta$ half-angle of the tip, $z$ piezo height, $d$ deflection, $k$ spring constant, $k_BT$ thermal energy.

## 1. Force curves

$$\delta = (z - z_0) - (d - d_0), \qquad F = k (d - d_0)$$

$(z_0, d_0)$ is the point of contact. The baseline is a straight line fitted on the non-contact part and subtracted.

Optional corrections: tilt (line fitted between two limits and subtracted) and smoothing (Savitzky-Golay filter [1]).

## 2. Point of contact

- **Detected**: last point of the curve below the baseline shifted by a threshold (residual noise of the baseline fit × a multiplier, 5 by default). The baseline can be fitted again just before the contact.
- **Fitted**: for a model $F = E g(\delta - \delta_0)$, $E$ has a closed form for each $\delta_0$:

$$E(\delta_0) = \frac{\sum F g}{\sum g^2}$$

The $\delta_0$ with the lowest error is searched on a grid, then refined (golden section). If the fit fails, the detected point is kept ("fixed point of contact").

## 3. Elasticity: contact models

**Paraboloid** (Hertz [2]), for $\delta \ll R$:

$$F = \frac{4}{3}\frac{E}{1-\nu^2}\sqrt{R} \delta^{3/2}$$

**Cone** (Sneddon [3]):

$$F = \frac{2}{\pi}\frac{E}{1-\nu^2}\tan\theta \delta^2$$

**Four-sided pyramid** (Barber & Billings [4], Rico et al. [5]), general form for $n$ sides with $n = 4$:

$$F = \frac{2\tan\theta}{n\sin(\pi/n)}\frac{E}{1-\nu^2} \delta^2 = \frac{1}{\sqrt{2}}\frac{E}{1-\nu^2}\tan\theta \delta^2$$

Within 6 % of the numerical solution of Bilodeau [6, 7], see also section 8.

**Flat punch** (Sneddon [3]):

$$F = 2R\frac{E}{1-\nu^2} \delta$$

**Sphere** (Sneddon [3]), with $a$ the contact radius:

$$\delta = \frac{a}{2}\ln\frac{R+a}{R-a}, \qquad F = \frac{E}{2(1-\nu^2)}\left[(a^2+R^2)\ln\frac{R+a}{R-a} - 2aR\right]$$

**Blunted pyramid** (Rico et al. [5], as written in [7]) and **blunted cone** (Briscoe et al. [8]; $c = 1$): paraboloid apex tangent to the faces at $b$.

$$c = \frac{n}{\pi}\sin\frac{\pi}{n}, \qquad b = \frac{cR}{\tan\theta}$$

For $a \le b$, Hertz. For $a > b$, $a$ is given by

$$\delta - \frac{ca}{\tan\theta}\left(\frac{\pi}{2} - \arcsin\frac{b}{a}\right) + \frac{a}{R}\left(\sqrt{a^2-b^2} - a\right) = 0$$

$$F = \frac{2E}{1-\nu^2}\left[\delta a - \frac{ca^2}{2\tan\theta}\left(\frac{\pi}{2} - \arcsin\frac{b}{a}\right) - \frac{a^3}{3R} + \sqrt{a^2-b^2}\left(\frac{cb}{2\tan\theta} + \frac{a^2-b^2}{3R}\right)\right]$$

**Needle** (Garcia & Garcia [9], eq. 18): flat punch capped with a cone, with $\delta_c = \pi R/(2\tan\theta)$:

$$F = \frac{2}{\pi}\frac{E}{1-\nu^2}\tan\theta \delta^2 \quad (\delta \le \delta_c), \qquad F = \frac{2}{\pi}\frac{E}{1-\nu^2}\tan\theta \delta_c^2 + 2R\frac{E}{1-\nu^2}(\delta-\delta_c) \quad (\delta > \delta_c)$$

**Nanowire** (Garcia & Garcia [9], eq. 17): flat punch capped with a hemisphere:

$$F = \frac{4}{3}\frac{E}{1-\nu^2}\sqrt{R} \delta^{3/2} \quad (\delta \le R), \qquad F = \frac{4}{3}\frac{E}{1-\nu^2}R^2 + 2R\frac{E}{1-\nu^2}(\delta-R) \quad (\delta > R)$$

**Thin shell** (Reissner [10]; Calò et al. [11]; Vorselen et al. [12]): point load on a shell of radius $R$ (the vesicle) and thickness $h$:

$$F = \frac{4Eh^2}{R\sqrt{3(1-\nu^2)}} \delta$$

## 4. Stiffness tomography (Slope)

Stiffness tomography (Roduit et al. [13], as used in Janel et al. [14]): slope $dF/d\delta$ of each segment of the force curve, by a linear fit, in N/m, without contact model.

## 5. Bottom effect correction

Garcia & Garcia [9], eqs. 12, 14 and 16 to 18, sample bonded to a rigid support of height $h$:

$$F = E g(\delta) \varphi(\chi)$$

- Paraboloid, $\chi = \sqrt{R\delta}/h$: $\varphi = 1 + 1.133\chi + 1.497\chi^2 + 1.469\chi^3 + 0.755\chi^4$
- Cone and pyramid, $\chi = \delta\tan\theta/h$: $\varphi = 1 + 0.721\chi + 0.650\chi^2 + 0.491\chi^3 + 0.225\chi^4$
- Flat punch, $\chi = R/h$: $\varphi = 1 + 1.133\chi + 1.283\chi^2 + 0.598\chi^3 - 0.291\chi^4$

**Nanowire** (eq. 17) and **needle** (eq. 18): the correction of the cap is applied up to the critical indentation $\delta_c$, where the contact radius reaches $R$, then the correction of the flat punch of radius $R$, with $\chi = R/h$ (constant).

Nanowire, paraboloid cap, $\delta_c = R$:

$$F = E g_{\text{par}}(\delta) \varphi_{\text{par}}(\delta) \quad (\delta \le R)$$

$$F = E g_{\text{par}}(R) \varphi_{\text{par}}(R) + 2R\frac{E}{1-\nu^2}(\delta - R) \varphi_{\text{punch}} \quad (\delta > R)$$

Needle, conical cap, $\delta_c = \pi R/(2\tan\theta)$:

$$F = E g_{\text{cone}}(\delta) \varphi_{\text{cone}}(\delta) \quad (\delta \le \delta_c)$$

$$F = E g_{\text{cone}}(\delta_c) \varphi_{\text{cone}}(\delta_c) + 2R\frac{E}{1-\nu^2}(\delta - \delta_c) \varphi_{\text{punch}} \quad (\delta > \delta_c)$$

$g_{\text{par}}$ and $g_{\text{cone}}$ are the force per unit modulus of the paraboloid and of the cone (section 3), and $\varphi_{\text{par}}$, $\varphi_{\text{cone}}$ and $\varphi_{\text{punch}}$ the polynomials above. The force is continuous at $\delta_c$, but not its slope. The influence of the support is small for thin tips: $\varphi_{\text{punch}} = 1.011$ for $R = 50$ nm and $h = 5$ µm.

The polynomials are expansions in the ratio of the contact radius to the height, and become unreliable when this ratio approaches 1.

The elasticity and the point of contact are fitted again with $\varphi$ in the model. The height is entered by hand, or is the topography above the plane of a substrate ROI.

## 6. Brush model

Sokolov et al. [15], Iyer et al. [16]. The cell body is fitted above a percentage of the maximum force (paraboloid or sphere), giving $E$ and the surface position $\delta_s$. The distance between the probe and the body is

$$h = \delta_s + \delta_{\text{body}}(F) - \delta$$

The brush force is fitted for $0.1 < h/L < 0.8$ (linear fit of $\ln F$ against $h$, iterated):

$$F = 100 k_BT R N^{3/2}L \exp\left(-\frac{2\pi h}{L}\right)$$

$L$ is the brush length and $N$ the grafting density. The prefactor is the one of Makarova et al. [17]; Guz et al. [18] use 50. No length is given for $R^2 < 0.8$.

## 7. Viscoelasticity: Ting fit

Power-law relaxation modulus (Efremov et al. [19, 20]), $t_0 = 1$ s, with fluidity $\beta$:

$$E(t) = E_0\left(\frac{t}{t_0}\right)^{-\beta}, \qquad \tan\delta_{\text{loss}} = \tan\frac{\pi\beta}{2}$$

Ting's solution [21], with $F = C E \delta^n$ for an elastic sample and $U = \delta^n$:

$$F(t) = C\int_0^{t} E(t-t')\frac{dU}{dt'} dt' \quad \text{(approach)}, \qquad F(t) = C\int_0^{t_1(t)} E(t-t')\frac{dU}{dt'} dt' \quad \text{(retraction)}$$

$$\int_{t_1}^{t} E(t-t')\frac{d\delta}{dt'} dt' = 0$$

At constant velocities $v_t$ and $v_r$, with $B$ the beta function and ${}_2F_1$ the hypergeometric function, $t_m$ the end of the approach:

$$F = C E_0 t_0^\beta v_t^n n B(n, 1-\beta) t^{n-\beta} \quad \text{(approach)}$$

$$F = C E_0 t_0^\beta v_t^n t_1^n t^{-\beta} {}_2F_1\left(\beta, n; n+1; \frac{t_1}{t}\right), \quad t_1 = t - \left(1+\frac{v_r}{v_t}\right)^{1/(1-\beta)}(t-t_m) \quad \text{(retraction)}$$

Otherwise the integrals are computed on the measured indentation. With the bottom effect correction, $U$ is multiplied by $\varphi$ (Tejedor & Garcia [22]). The fit routine follows TingFit of PyFMRheo [23].

## 8. Microrheology

Alcaraz et al. [24], PyFMRheo [23]. With $H(f) = F(f)/\delta(f)$, by FFT or sine fit, optionally corrected with a piezo characterization, and the contact model linearized around the working indentation $\delta_0$:

$$G^*(f) = G'(f) + iG''(f) = \frac{H(f) - i 2\pi f b(0)}{2(1+\nu) g'(\delta_0)}$$

For a paraboloid:

$$G^*(f) = \frac{1-\nu}{4\sqrt{R\delta_0}}\left[H(f) - i 2\pi f b(0)\right]$$

Viscous drag (Alcaraz et al. [25]), measured above the surface at several heights and extrapolated to the surface:

$$b(h) = \frac{\mathrm{Im}(F/\delta)}{2\pi f}, \qquad b(h) = \frac{6\pi\eta a_{\text{eff}}^2}{h + h_{\text{eff}}}, \qquad b(0) = \frac{6\pi\eta a_{\text{eff}}^2}{h_{\text{eff}}}$$

For the four-sided pyramid, Alcaraz et al. [24] use $F = \frac{3E\tan\theta}{4(1-\nu^2)}\delta^2$ (their eq. 1) and $G^* = \frac{1-\nu}{3\delta_0\tan\theta}\left[H - i2\pi f b(0)\right]$ (their eq. 6). pyAF uses the coefficient of section 3, $1/\sqrt{2}$ instead of $3/4$, which gives moduli about 6 % higher.

## 9. Adhesion

Both are measured from the line through the two jumps of contact of the retract curve:

- **Work of detachment**: area between this line and the curve (trapezes).
- **Detachment force**: largest distance between this line and the curve.

## 10. Unbinding events

- **Kernel**: convolution of the curve with a step kernel, normalized on the noise, with an optional adaptive threshold (a Savitzky-Golay smoothing of the convolution is subtracted).
- **MSF**: in a sliding window, comparison of a step fit and of a single line: the difference of their residuals, multiplied by the height of the step and normalized.
- **Force of an event**: $|y_1 - y_2|$, the fit through the step crossing the linear fits on each side.
- **Loading rate**: $r = \frac{dF}{dz} v$, slope over the rise of the force after the event, up to a fraction of the force of the event.

## 11. Dynamic force spectroscopy

See Mesbah et al. [26]. All models are fitted on median forces.

**Bell-Evans** (Evans & Ritchie [27]):

$$F = \frac{k_BT}{x}\ln\left(1 + \ln 2 \frac{r x}{k_BT k_{\text{off}}}\right)$$

**Dudko-Hummer-Szabo** [28], $\nu = 1/2$ (cusp) or $2/3$ (linear-cubic):

$$F = \frac{\Delta G}{\nu x}\left(1 - \left[1 - \frac{k_BT}{\Delta G}\ln\left(1 + \ln 2 \frac{r x}{k_BT k_{\text{off}}}\right)\right]^{\nu}\right)$$

**Friddle-Noy-De Yoreo** [29], with rebinding:

$$F = f_{eq} + \frac{k_BT}{x}\ln\left(1 + \ln 2 \frac{r x}{k(f_{eq}) k_BT}\right), \quad k(f_{eq}) = k_{\text{off}} e^{f_{eq}x/k_BT}, \quad f_{eq} = \sqrt{2k_c\Delta G_{bu}}$$

The original formulas give the most probable force (Bell-Evans) or the mean force (Dudko-Hummer-Szabo, Friddle-Noy-De Yoreo). pyAF fits median forces, less affected by multiple bonds and nonspecific adhesions, so it uses the corresponding expression, with $\ln 2$.

## 12. Breakthrough of lipid bilayers

See Gumí-Audenis & Giannotti [30]. A breakthrough is a jump of the separation $z - d$ above its local trend, larger than a minimum thickness and $n_\sigma$ times the noise, with a drop of the force. $F_b$ is the force just before the jump, and the thickness is the size of the jump.

## 13. Statistics

- **SuperPlots**: one value per file, the median of its curves (the mean for the events per curve), and statistics on the replicates.
- **Tests on the means of the replicates**: t-test, paired or not, Tukey HSD for all the pairs, Holm correction, Fieller interval for the fold change.
- **Mixed-effects model on the curves**: value = mean of the condition + offset of the replicate + offset of the file + residual, by REML (statsmodels).
- **Mode**: maximum of a Gaussian kernel density estimate (Scott's rule).

## References

1. Savitzky A, Golay MJE. Smoothing and differentiation of data by simplified least squares procedures. *Anal. Chem.* 36:1627–1639 (1964). [doi:10.1021/ac60214a047](https://doi.org/10.1021/ac60214a047)
2. Hertz H. Ueber die Berührung fester elastischer Körper. *J. Reine Angew. Math.* 92:156–171 (1882). [doi:10.1515/crll.1882.92.156](https://doi.org/10.1515/crll.1882.92.156)
3. Sneddon IN. The relation between load and penetration in the axisymmetric Boussinesq problem for a punch of arbitrary profile. *Int. J. Eng. Sci.* 3:47–57 (1965). [doi:10.1016/0020-7225(65)90019-4](https://doi.org/10.1016/0020-7225(65)90019-4)
4. Barber JR, Billings DA. An approximate solution for the contact area and elastic compliance of a smooth punch of arbitrary shape. *Int. J. Mech. Sci.* 32:991–997 (1990). [doi:10.1016/0020-7403(90)90003-2](https://doi.org/10.1016/0020-7403(90)90003-2)
5. Rico F, Roca-Cusachs P, Gavara N, et al. Probing mechanical properties of living cells by atomic force microscopy with blunted pyramidal cantilever tips. *Phys. Rev. E* 72:021914 (2005). [doi:10.1103/PhysRevE.72.021914](https://doi.org/10.1103/PhysRevE.72.021914)
6. Bilodeau GG. Regular pyramid punch problem. *J. Appl. Mech.* 59:519–523 (1992). [doi:10.1115/1.2893754](https://doi.org/10.1115/1.2893754)
7. Lekka M, Navajas D, Radmacher M, Podestà A (eds). *Mechanics of Cells and Tissues in Diseases, Volume 1: Biomedical Methods*, section 2.2, table 2.2.1. De Gruyter (2023). [doi:10.1515/9783110640632](https://doi.org/10.1515/9783110640632)
8. Briscoe BJ, Sebastian KS, Adams MJ. The effect of indenter geometry on the elastic response to indentation. *J. Phys. D: Appl. Phys.* 27:1156–1162 (1994). [doi:10.1088/0022-3727/27/6/013](https://doi.org/10.1088/0022-3727/27/6/013)
9. Garcia PD, Garcia R. Determination of the elastic moduli of a single cell cultured on a rigid support by force microscopy. *Biophys. J.* 114:2923–2932 (2018). [doi:10.1016/j.bpj.2018.05.012](https://doi.org/10.1016/j.bpj.2018.05.012)
10. Reissner E. Stresses and small displacements of shallow spherical shells. *J. Math. Phys.* 25:80–85 and 279–300 (1946). [doi:10.1002/sapm194625180](https://doi.org/10.1002/sapm194625180), [doi:10.1002/sapm1946251279](https://doi.org/10.1002/sapm1946251279)
11. Calò A, Reguera D, Oncins G, et al. Force measurements on natural membrane nanovesicles reveal a composition-independent, high Young's modulus. *Nanoscale* 6:2275 (2014). [doi:10.1039/C3NR05107B](https://doi.org/10.1039/C3NR05107B)
12. Vorselen D, Piontek MC, Roos WH, et al. Mechanical characterization of liposomes and extracellular vesicles, a protocol. *Front. Mol. Biosci.* 7:139 (2020). [doi:10.3389/fmolb.2020.00139](https://doi.org/10.3389/fmolb.2020.00139)
13. Roduit C, Sekatski S, Dietler G, et al. Stiffness tomography by atomic force microscopy. *Biophys. J.* 97:674–677 (2009). [doi:10.1016/j.bpj.2009.05.010](https://doi.org/10.1016/j.bpj.2009.05.010)
14. Janel S, Popoff M, Barois N, et al. Stiffness tomography of eukaryotic intracellular compartments by atomic force microscopy. *Nanoscale* 11:10320–10328 (2019). [doi:10.1039/C8NR08955H](https://doi.org/10.1039/C8NR08955H)
15. Sokolov I, Iyer S, Subba-Rao V, et al. Detection of surface brush on biological cells in vitro with atomic force microscopy. *Appl. Phys. Lett.* 91:023902 (2007). [doi:10.1063/1.2757104](https://doi.org/10.1063/1.2757104)
16. Iyer S, Gaikwad RM, Subba-Rao V, et al. Atomic force microscopy detects differences in the surface brush of normal and cancerous cells. *Nat. Nanotechnol.* 4:389–393 (2009). [doi:10.1038/nnano.2009.77](https://doi.org/10.1038/nnano.2009.77)
17. Makarova N, Lekka M, Gnanachandran K, et al. Mechanical way to study molecular structure of pericellular layer. *ACS Appl. Mater. Interfaces* 15:35962–35972 (2023). [doi:10.1021/acsami.3c06341](https://doi.org/10.1021/acsami.3c06341)
18. Guz N, Dokukin M, Kalaparthi V, et al. If cell mechanics can be described by elastic modulus: study of different models and probes used in indentation experiments. *Biophys. J.* 107:564–575 (2014). [doi:10.1016/j.bpj.2014.06.033](https://doi.org/10.1016/j.bpj.2014.06.033)
19. Efremov YM, Wang WH, Hardy SD, et al. Measuring nanoscale viscoelastic parameters of cells directly from AFM force-displacement curves. *Sci. Rep.* 7:1541 (2017). [doi:10.1038/s41598-017-01784-3](https://doi.org/10.1038/s41598-017-01784-3)
20. Efremov YM, Okajima T, Raman A. Measuring viscoelasticity of soft biological samples using atomic force microscopy. *Soft Matter* 16:64–81 (2020). [doi:10.1039/C9SM01020C](https://doi.org/10.1039/C9SM01020C)
21. Ting TCT. The contact stresses between a rigid indenter and a viscoelastic half-space. *J. Appl. Mech.* 33:845–854 (1966). [doi:10.1115/1.3625192](https://doi.org/10.1115/1.3625192)
22. Tejedor JR, Garcia R. High-throughput nanorheology of living cells powered by supervised machine learning. *Adv. Intell. Syst.* 7:2400867 (2025). [doi:10.1002/aisy.202400867](https://doi.org/10.1002/aisy.202400867)
23. PyFMRheo (PyFMLab), DyNaMo, INSERM. [github.com/DyNaMo-INSERM/PyFMlab_DyNaMo](https://github.com/DyNaMo-INSERM/PyFMlab_DyNaMo)
24. Alcaraz J, Buscemi L, Grabulosa M, et al. Microrheology of human lung epithelial cells measured by atomic force microscopy. *Biophys. J.* 84:2071–2079 (2003). [doi:10.1016/S0006-3495(03)75014-0](https://doi.org/10.1016/S0006-3495(03)75014-0)
25. Alcaraz J, Buscemi L, Puig-de-Morales M, et al. Correction of microrheological measurements of soft samples with atomic force microscopy for the hydrodynamic drag on the cantilever. *Langmuir* 18:716–721 (2002). [doi:10.1021/la0110850](https://doi.org/10.1021/la0110850)
26. Mesbah I, Sumbul F, Valotteau C, Rico F. Single-molecule force spectroscopy: experiments, analysis, and simulations. In Santos NC, Carvalho FA (eds), *Atomic Force Microscopy: Methods and Protocols*, 2nd ed., *Methods Mol. Biol.* vol. 3077, chap. 5, pp. 105–136 (2026). [doi:10.1007/978-1-0716-5555-9_5](https://doi.org/10.1007/978-1-0716-5555-9_5)
27. Evans E, Ritchie K. Dynamic strength of molecular adhesion bonds. *Biophys. J.* 72:1541–1555 (1997). [doi:10.1016/S0006-3495(97)78802-7](https://doi.org/10.1016/S0006-3495(97)78802-7)
28. Dudko OK, Hummer G, Szabo A. Intrinsic rates and activation free energies from single-molecule pulling experiments. *Phys. Rev. Lett.* 96:108101 (2006). [doi:10.1103/PhysRevLett.96.108101](https://doi.org/10.1103/PhysRevLett.96.108101)
29. Friddle RW, Noy A, De Yoreo JJ. Interpreting the widespread nonlinear force spectra of intermolecular bonds. *Proc. Natl. Acad. Sci. USA* 109:13573–13578 (2012). [doi:10.1073/pnas.1202946109](https://doi.org/10.1073/pnas.1202946109)
30. Gumí-Audenis B, Giannotti MI. Structural and mechanical characterization of supported model membranes by AFM. In *Biomimetic Lipid Membranes: Fundamentals, Applications, and Commercialization*, pp. 1–27. Springer (2019). [doi:10.1007/978-3-030-11596-8_1](https://doi.org/10.1007/978-3-030-11596-8_1)
