# Exact counterexamples and stable oscillations in a two-host dormancy model

## Abstract

We study the six-dimensional deterministic host–virus system introduced by Blath and Tóbiás (2026). Their Conjecture 4.1(C)(iii) predicts stable positive coexistence when both resident endemic subsystems are locally asymptotically stable and each host can invade the other. We give an exact one-parameter family satisfying these hypotheses whose unique positive equilibrium has a two-dimensional unstable subspace. Every nonnegative equilibrium is hyperbolic and unstable, whereas all positive solutions are global and bounded. Consequently, the positive initial conditions that converge to any equilibrium form a set of Lebesgue measure zero. We then fix all parameters except the dormancy probability and certify a supercritical Hopf bifurcation, producing orbitally asymptotically stable positive periodic orbits. The proof combines explicit elimination, rational Hurwitz certificates, a stable-manifold argument, and an outward integer-interval calculation of the Hopf crossing and first Lyapunov coefficient. The portable supplement reconstructs the Jacobian from the published equations, verifies polynomial identities exactly, and certifies the normal-form signs. The result concerns deterministic long-time behaviour; it does not contradict the source's proved stochastic invasion results.

## 1. The question and the model

Contact-mediated dormancy changes both the immediate effect of viral contact and the subsequent supply of susceptible hosts. The model of Blath and Tóbiás (2026) contains a host without dormancy and a competing host that may become dormant following viral contact. The source develops stochastic early-invasion results and formulates separate conjectures about the deterministic system's eventual behaviour. Conjecture 4.1(C)(iii) links stable resident endemic states and strict mutual invasion to positive coexistence stability and attraction. Here we preserve those hypotheses, disprove the proposed stability conclusion, and identify a local family of attracting periodic coexistence states. Invasion growth and coexistence stability concern different linearised systems.

Let $a,i$ be the active and infected densities of host 1; $b,d,j$ the active, dormant and infected densities of host 2; and $V$ the free viral density. Put $H=a+i+b+d+j$. In the state order $X=(a,i,b,d,j,V)$, Eq. (1.1) of Blath and Tóbiás (2026) is

$$
\begin{aligned}
\dot a&=a(\lambda_1-\mu_1-CH-DV)+ri,\\
\dot i&=DVa-(r+v)i,\\
\dot b&=b(\lambda_2-\mu_1-CH-DV)+rj+\sigma d,\\
\dot d&=qDVb-(\kappa\mu_1+\sigma)d,\\
\dot j&=(1-q)DVb-(r+v)j,\\
\dot V&=-DV[a+(1-q)b]+mv(i+j)-\mu_3V.
\end{aligned} \tag{1.1}
$$

Here $q$ is the probability of entering dormancy. One preceding prose description in the source reverses $q$ and $1-q$; the equations determine our convention. The admissible parameters satisfy $0<\mu_1<\lambda_2<\lambda_1$, $0<q<1$, positive $C,D,r,v,\sigma,\mu_3$, nonnegative $\kappa$, and integer burst size $m$. Our examples have $\kappa>0$ and $m=3$.

A matrix is *Hurwitz stable* when every eigenvalue has strictly negative real part. Resident stability means stability within the subsystem containing one host and virus; it does not mean stability against addition of the missing host. Earlier work treats contact-mediated dormancy and one-host virus–dormancy systems (Gulbudak & Weitz, 2016; Blath & Tóbiás, 2023). The contribution here is the explicit counterexamples and certified single-parameter bifurcation with the source hypotheses preserved. The general instability and Hopf methods are classical.

## 2. An exact family of counterexamples

For every $C\geq\tfrac12$, choose

$$
\begin{aligned}
&\mu_1=D=v=1,\qquad r=\frac12,\qquad q=\frac{49}{50},\\
&\kappa=\frac1{20},\qquad \sigma=\frac9{20},\qquad m=3,\qquad\mu_3=\frac75,\\
&\lambda_1=\frac{13+665C}{3},\qquad
\lambda_2=\frac{467}{300}+\frac{665C}{3}.
\end{aligned} \tag{2.1}
$$

**Theorem 2.1 (exact family).** For (1.1) with (2.1), both resident endemic equilibria exist and are Hurwitz stable in their resident subsystems. Both invasion growth rates are strictly positive. The unique strictly positive equilibrium is

$$
X_*=\left(1,\frac{10}{3},20,196,\frac43,5\right). \tag{2.2}
$$

Its Jacobian has exactly two eigenvalues in the open right half-plane and four in the open left half-plane. The two unstable eigenvalues are nonreal. There are exactly six nonnegative equilibria, all hyperbolic and unstable in the full system. Every strictly positive solution exists globally and is bounded. The strictly positive initial states that converge to any equilibrium form a set of six-dimensional Lebesgue measure zero.

The strict hypotheses and the failure of stable positive coexistence persist on an open set of nearby admissible parameters around each family member, retaining $m=3$.

This refutes the deterministic assertion of Conjecture 4.1(C)(iii). The almost-everywhere conclusion concerns convergence to equilibria. It does not identify the eventual attractor for each $C$ or assert uniform persistence of every coordinate.

### 2.1. Equilibrium elimination and the source hypotheses

The construction can be read backwards from a positive equilibrium. Fixing $a=1$, $b=20$ and $V=5$ determines the infected and dormant populations through their linear balance equations. The virus balance then selects $\mu_3=7/5$, and the two active-host balances determine the affine birth rates in (2.1). Their common slope $665/3$ is the total equilibrium host density. Consequently varying $C$ while adjusting both birth rates preserves the same positive equilibrium; it does not represent changing competition alone. The remaining task is to verify resident stability, mutual invasion and coexistence instability simultaneously, which the exact sign calculations below do for the entire ray. This is a reproducible algebraic construction, not an assertion about how the original parameter search was conducted.

The infected and dormant balances at a positive equilibrium give

$$
i=\frac23Va,\qquad d=\frac{49}{25}Vb,\qquad j=\frac1{75}Vb.
$$

Dividing the two active-host balances by their positive active densities and subtracting yields

$$
\lambda_1-\lambda_2=\frac{833}{1500}V,
\qquad \lambda_1-\lambda_2=\frac{833}{300}.
$$

Thus $V=5$. Virus balance reduces to $a+b/50=7/5$, and the first active-host balance gives $H=665/3$. Hence

$$
\begin{aligned}
&i=\frac{10}{3}a,\qquad d=\frac{49}{5}b,\qquad j=\frac1{15}b,\\
&\frac{13}{3}a+\frac{163}{15}b=\frac{665}{3},\qquad
a+\frac{b}{50}=\frac75.
\end{aligned}
$$

This nonsingular linear system gives $a=1,b=20$, proving existence, uniqueness and positivity of (2.2). The source's nondegeneracy condition holds because

$$
r\kappa\mu_1=\frac1{40}\ne\frac9{20}=v\sigma.
$$

All basic parameter inequalities hold, and $mv=3>r+v=3/2$. The active resident densities are strictly below the virus-free carrying densities:

$$
\frac75<\frac{\lambda_1-1}{C},\qquad
70<\frac{\lambda_2-1}{C}.
$$

Indeed, the two quantities on the right are $665/3+10/(3C)$ and $665/3+167/(300C)$. These inequalities explicitly check endemic resident existence despite an inconsistency in the source's numbering of the corresponding conditions.

### 2.2. Resident equilibria and strict mutual invasion

The two resident endemic states are

$$
\begin{aligned}
E_1&=\left(\frac75,\frac{14V_1}{15},0,0,0,V_1\right),\\
E_2&=\left(0,0,70,\frac{686V_2}{5},\frac{14V_2}{15},V_2\right),
\end{aligned} \tag{2.3}
$$

where

$$
V_1=\frac{25+1652C}{5+7C},\qquad
V_2=\frac{5(167+45500C)}{167+207200C}. \tag{2.4}
$$

Substitution verifies every equilibrium equation. Appendix A gives the resident cubic and quartic and the exact Hurwitz certificates establishing resident stability for every $C\geq1/2$.

We also determine the exact unstable dimensions of the missing-host blocks. At $E_1$, the invader coordinates $(b,d,j)$ have matrix

$$
M_2=
\begin{pmatrix}
\lambda_2-1-CH_1-V_1&9/20&1/2\\
(49/50)V_1&-1/2&0\\
V_1/50&0&-3/2
\end{pmatrix}, \tag{2.5}
$$

where $H_1$ is the resident host total. Positive diagonal similarity makes this arrowhead matrix symmetric, since both pairs of opposite off-diagonal entries are positive. A congruence eliminating its two negative diagonal entries leaves scalar Schur complement

$$
G_2=\frac{833}{1500}(V_1-5).
$$

Sylvester's law of inertia therefore gives one positive and two negative eigenvalues when $G_2>0$. At $E_2$, the analogous block on $(a,i)$ is

$$
M_1=\begin{pmatrix}
\lambda_1-1-CH_2-V_2&1/2\\
V_2&-3/2
\end{pmatrix}. \tag{2.6}
$$

It is likewise similar to a real symmetric matrix, with scalar Schur complement $G_1=(833/1500)(5-V_2)$. Both invasion tests are strict:

$$
V_1-5=\frac{1617C}{5+7C}>0,\qquad
5-V_2=\frac{808500C}{167+207200C}>0. \tag{2.7}
$$

Each resident endemic equilibrium consequently has exactly one unstable eigenvalue in the full system. All of its other eigenvalues have negative real part.

### 2.3. Exact instability at coexistence

Direct differentiation gives $p_C(s)=\det(sI-J(X_*))$:

$$
\begin{aligned}
18000p_C(s)={}&18000s^6+(223380+378000C)s^5\\
&+(782884+4951380C)s^4+(746367+18467484C)s^3\\
&+(525350+19832547C)s^2+(673869+11876130C)s\\
&+12122649C.
\end{aligned} \tag{2.8}
$$

For descending coefficients $a_0,\ldots,a_n$, define the leading Hurwitz determinants by

$$
\Delta_k=\det[a_{2j-i}]_{i,j=1}^k,\qquad
a_h=0\quad\hbox{if }h\notin\{0,\ldots,n\}.
$$

For (2.8), their signs are

$$
(\operatorname{sgn}\Delta_1,\ldots,\operatorname{sgn}\Delta_6)
=(+,+,+,+,-,-). \tag{2.9}
$$

An exact coefficient certificate covers the entire unbounded ray. Substitute $t=2C-1\geq0$ and clear positive denominators. Every coefficient of each of the first four determinant numerators is strictly positive; every coefficient of the fifth and sixth is strictly negative. The degrees are $1,2,3,4,5,6$. The supplement supplies every coefficient and verifies the determinant expansion with both a symbolic rational implementation and a package-free integer polynomial implementation.

The Routh first-column signs are $(+,+,+,+,+,-,+)$, with two sign changes and no zero pivot. The Routh criterion gives exactly two right-half-plane roots and no imaginary-axis roots. All coefficients in (2.8) are positive, so no positive real number is a root. The unstable roots must form a nonreal conjugate pair.

At the concrete value $C=1/2$, the Jacobian is

$$
J_*=
\begin{pmatrix}
-13/6&0&-1/2&-1/2&-1/2&-1\\
5&-3/2&0&0&0&1\\
-10&-10&-4333/300&-191/20&-19/2&-20\\
0&0&49/10&-1/2&0&98/5\\
0&0&1/10&0&-3/2&2/5\\
-5&3&-1/10&0&3&-14/5
\end{pmatrix}. \tag{2.10}
$$

Its trace is $-2291/100=-22.91$. The unstable roots are approximately $0.02111770\pm0.84217481i$; the proof uses exact signs.

### 2.4. All equilibria and boundedness

If $V=0$, equilibrium forces $i=d=j=0$. Since $\lambda_1\ne\lambda_2$, the two active hosts cannot both be present. This gives the origin and the two virus-free single-host equilibria. If $V>0$, at least one active host is positive; the three possibilities yield exactly $E_1,E_2,X_*$. A missing active host has zero infected density, and a missing host 2 has zero dormant density. Thus no additional face equilibrium is omitted.

Write $K_h=(\lambda_h-1)/C$. At the origin the two active-host eigenvalues $\lambda_h-1$ are positive and the other four eigenvalues are negative. At the virus-free host-1 state, the active resident direction, the missing host-2 direction, and its dormant and infected directions are stable. The infection–virus block is

$$
\begin{pmatrix}-3/2&K_1\\3&-K_1-7/5\end{pmatrix},
$$

with determinant $21/10-3K_1/2<0$. It has one positive and one negative eigenvalue. At the virus-free host-2 state, the missing active host 1 grows at rate $\lambda_1-\lambda_2>0$. The host-2 infection–virus block has determinant $21/10-3(1-q)K_2/2<0$, since $K_2>70$, and the remaining directions are stable. The six full-system unstable dimensions are therefore:

| Equilibrium | Unstable dimension |
| --- | ---: |
| Origin | 2 |
| Virus-free host 1 | 1 |
| Virus-free host 2 | 2 |
| Endemic host 1 | 1 |
| Endemic host 2 | 1 |
| Positive coexistence | 2 |

All six equilibria are hyperbolic. Local resident stability and full-system instability have been checked in their distinct spaces.

The field is quasipositive, and a positive solution remains strictly positive on its interval of existence. For boundedness take $L=4H+V$, $A=\lambda_1-1$ and $\delta=1/40$. Summation gives

$$
H'=a(\lambda_1-1-CH)+b(\lambda_2-1-CH)-(i+j)-\frac1{20}d.
$$

With $s=a+b\leq H$,

$$
\begin{aligned}
L'+\delta L
&\leq4(A+\delta)s-4CHs\\
&\leq4(A+\delta)s-4Cs^2
\leq\frac{(A+\delta)^2}{C}.
\end{aligned} \tag{2.11}
$$

The omitted coefficients of $i,j,d,V$ on the first line are nonpositive, as is the viral adsorption term. Scalar comparison bounds $L$ for all positive times; polynomial local existence then extends each positive solution globally.

### 2.5. Almost-everywhere nonconvergence to equilibria

For each equilibrium, the stable-manifold theorem gives a local stable manifold of dimension less than six. A convergent trajectory eventually belongs to that local stable manifold. Its basin is therefore contained in the union, over integer times, of finite-time inverse images of that manifold. On the open positive cone each finite-time flow is a local smooth diffeomorphism onto its image, so every such inverse image has Lebesgue measure zero. The countable union over times and finite union over equilibria still have measure zero.

Every property needed for the counterexample is strict: the face equilibria are positive in their respective faces, the resident and invasion conditions are strict, and the coexistence equilibrium is hyperbolic with an unstable pair. The equilibrium elimination is nonsingular. Continuity preserves these properties on an open admissible parameter neighbourhood of each family member with $m=3$ fixed. This completes Theorem 2.1. $\square$

## 3. Attracting periodic coexistence with only dormancy varied

Now fix $C=1/2$ and all parameters in (2.1) except $q$. In particular, both birth rates remain fixed:

$$
\lambda_1=\frac{691}{6},\qquad \lambda_2=\frac{11239}{100}.
$$

The positive equilibrium branch is

$$
\begin{aligned}
V(q)&=\frac{49}{10q},&
b(q)&=\frac{17020q-833}{15q(5q+49)},\\
a(q)&=\frac75-(1-q)b(q),&
i(q)&=\frac23V(q)a(q),\\
d(q)&=2qV(q)b(q),&
j(q)&=\frac23(1-q)V(q)b(q).
\end{aligned} \tag{3.1}
$$

**Theorem 3.1 (q-only Hopf bifurcation).** Let $P$ be the primitive degree-17 polynomial defined below. It has a unique root $q_H$ between the exact rational endpoints

$$
\begin{aligned}
q_-&=0.9769161486447415813763449883881468671926,\\
q_+&=0.9769161486447415813763449883881468671927.
\end{aligned} \tag{3.2}
$$

The branch (3.1) undergoes a nondegenerate supercritical Hopf bifurcation at $q_H$. For $q>q_H$ sufficiently close to $q_H$, there are orbitally asymptotically stable, strictly positive periodic orbits. Both resident endemic subsystems are stable and both invasion growth rates are positive in a neighbourhood of this bifurcation.

Orbital stability permits a phase shift along the periodic orbit. The neighbourhood is existential: the theorem does not certify a numerical radius, a periodic orbit at $q=0.98$, or global continuation of the periodic branch.

![Numerical spectral abscissa at coexistence: (a) the family varying C; (b) the branch varying only q. The marker is the rigorously isolated Hopf point; shading denotes equilibrium instability. Attracting periodic orbits are proved only sufficiently near the Hopf point.](../figures/dormancy_spectral_abscissa.png)

### 3.1. Exact isolation of the imaginary pair

Here the continuation differs from Section 2: both birth rates and $C=1/2$ remain fixed, and only $q$ varies. The explicit equilibrium branch permits substitution into the characteristic coefficients before taking a resultant or Hurwitz determinant. Isolating a zero of that determinant is only a candidate-location step; the complementary-root, transversality and nonlinear checks below are what turn it into a Hopf certificate.

Let $A(q)$ be the Jacobian at (3.1), with characteristic polynomial

$$
\det(zI-A(q))=z^6+\alpha_1z^5+\alpha_2z^4+\alpha_3z^3
+\alpha_4z^2+\alpha_5z+\alpha_6.
$$

Appendix B gives all six rational coefficient functions. Let $D_k$ be the associated Hurwitz determinants. Exact symbolic calculation proves

$$
D_5(q)=\frac{-343P(q)}{295245\,10^{18}q^9(5q+49)^8}, \tag{3.3}
$$

where $P$ is primitive with positive leading coefficient. Equation (3.3) and Appendix B define $P$ exactly; its eighteen integer coefficients are also supplied separately in ascending order.

Integer evaluation gives $P(q_-)<0<P(q_+)$. If $P(q)=\sum_{k=0}^{17}P_kq^k$, then on $[0,1]$,

$$
|P''(q)|\leq M_2:=\sum_{k=0}^{17}k(k-1)|P_k|.
$$

For $q_0=(q_-+q_+)/2$, exact rational arithmetic verifies

$$
P'(q_0)>M_2(q_+-q_-)/2>0.
$$

Hence $P'>0$ throughout the interval. The intermediate value theorem gives a unique root there.

Put

$$
\begin{aligned}
T&=\alpha_1\alpha_2\alpha_3-\alpha_3^2
+\alpha_1\alpha_5-\alpha_1^2\alpha_4,\\
x&=\frac{\alpha_1\alpha_2\alpha_5-\alpha_3\alpha_5-\alpha_1^2\alpha_6}{T}.
\end{aligned} \tag{3.4}
$$

Here $T=D_3$. Exact polynomial division verifies that the numerators of

$$
\alpha_1x^2-\alpha_3x+\alpha_5,\qquad
-x^3+\alpha_2x^2-\alpha_4x+\alpha_6
$$

are divisible by $P$. Their denominators are nonzero on the certified interval. Thus both expressions vanish at $q_H$. With $\omega=\sqrt{x}>0$, the characteristic polynomial factors at $q_H$ as

$$
(z^2+x)(z^4+\beta_1z^3+\beta_2z^2+\beta_3z+\beta_4), \tag{3.5}
$$

where

$$
\beta_1=\alpha_1,\quad \beta_2=\alpha_2-x,\quad
\beta_3=\alpha_3-\alpha_1x,\quad \beta_4=\alpha_6/x.
$$

The supplement verifies this factorisation coefficient by coefficient modulo $P$. A vanishing Hurwitz determinant alone would not suffice to establish a simple imaginary pair or stability of the remaining eigenvalues; (3.4)–(3.5) provide that additional structure.

### 3.2. Outward interval verification

The interval computation evaluates the equilibrium, $T$, $x$ and the quartic over $[q_-,q_+]$. Every equilibrium coordinate, $T$, $x$, and each quartic Hurwitz determinant is strictly positive. Coarse enclosures contained in the full certificate are:

| Quantity | Certified enclosure |
| --- | --- |
| $\beta_1$ | $(22.8908,22.8910)$ |
| $\beta_2$ | $(180.1628,180.1630)$ |
| $\beta_3$ | $(538.9358,538.9361)$ |
| $\beta_4$ | $(464.8550,464.8554)$ |
| Quartic $\Delta_2$ | $(3585.14,3585.16)$ |
| Quartic $\Delta_3$ | $(1688584,1688586)$ |
| Quartic $\Delta_4$ | $(784947289,784947292)$ |

Consequently, $\pm i\omega$ are simple and the other four eigenvalues lie in the open left half-plane. The certified frequency, total crossing speed and first Lyapunov coefficient satisfy

$$
0.824358242148467<\omega_H<0.824358242148468, \tag{3.6}
$$

$$
7.005341843234863<
\left.\frac{d}{dq}\operatorname{Re}\lambda(q)\right|_{q_H}
<7.005341843234864, \tag{3.7}
$$

$$
-0.024837104528915<l_1<-0.024837104528914. \tag{3.8}
$$

The derivative differentiates $A(q)$ along the moving equilibrium branch; it does not hold the equilibrium state fixed.

Choose right and transpose-left eigenvectors satisfying

$$
Au=i\omega u,\qquad \ell^TA=i\omega\ell^T,
\qquad u_6=1,\qquad \ell^Tu=1.
$$

Then the crossing speed is $\operatorname{Re}(\ell^TA'u)$. This convention is equivalent to the usual conjugate-left convention after conjugating the left eigenvector. Since the field is quadratic in the state, $D^3F=0$. With $B=D^2F$ from Appendix C, our normalisation is

$$
\begin{aligned}
l_1=\frac1{2\omega}\operatorname{Re}\Bigl[\ell^T\bigl\{
&-2B\bigl(u,A^{-1}B(u,\bar u)\bigr)\\
&+B\bigl(\bar u,(2i\omega I-A)^{-1}B(u,u)\bigr)
\bigr\}\Bigr].
\end{aligned} \tag{3.9}
$$

This is the standard Hopf normal-form expression in the stated normalisation (Kuznetsov, 1998). Its magnitude depends on the amplitude normalisation of $u$, whereas its sign does not.

The implementation represents a real interval by two integers on the grid $10^{-60}\mathbb Z$. Addition and subtraction use exact endpoint operations; multiplication and division round outwards with integer floor and ceiling; square roots use integer square roots. Complex intervals are rectangles. Every division excludes zero from its denominator enclosure. Integer midpoint scores only order candidate pivots in Gaussian elimination. A pivot is accepted only after interval verification that it is nonzero. No floating-point eigenvalue or midpoint supplies a proof bound.

The first five equations of $(A-i\omega I)u=0$ determine $u$ with $u_6=1$. The certified nonzero principal minor and exact vanishing of the full determinant imply the omitted equation at $q_H$. The transpose construction gives a left eigenvector; its pairing with $u$ is separated from zero before normalisation. Further certified linear solves evaluate (3.9). Interval automatic differentiation of the branch and the Jacobian gives $A'$, while the symbolic supplement records its full exact formula.

Fresh execution yields enclosure widths below $10^{-32}$ for the crossing and $10^{-31}$ for $l_1$, much smaller than the displayed sign margins. The full endpoints and every elimination pivot are recorded. Eigenvector residual intervals contain zero as diagnostics; exact eigenvector existence follows from the algebraic rank argument, not from residual containment.

### 3.3. Completion of the Hopf argument

The field is smooth in state and parameter near the interior equilibrium. The simple imaginary pair, stable quartic, positive crossing and negative nonzero $l_1$ meet the classical Hopf theorem. It yields a supercritical periodic branch for $q>q_H$. The four complementary directions and the Hopf radial direction are attracting, giving orbital asymptotic stability. The cycles shrink to a strictly positive equilibrium as $q\downarrow q_H$, so all their coordinates are positive when they are sufficiently small.

The interval supplement also checks both resident Hurwitz criteria and both invasion Schur complements over $[q_-,q_+]$. The latter lie in $(2.0265,2.0267)$ and $(52.6470,52.6471)$. These quantities determine the signs of the dominant invasion eigenvalues by the arrowhead argument in Section 2.2; they are not themselves eigenvalues or exponential growth rates. Continuity preserves these strict properties in an open neighbourhood of $q_H$. This proves Theorem 3.1. $\square$

## 4. Interpretation, scope and reproducibility

The two parameter paths answer different questions. In Theorem 2.1, competition strength and both birth rates vary together to keep the coexistence state fixed. In Theorem 3.1, the birth rates and competition strength remain fixed while dormancy changes and the equilibrium moves. These rational parameters constitute a mathematical test of a universal model conjecture; they have not been fitted to microbial measurements. Resident invasion examines the growth of a rare missing host in a fixed resident environment. At coexistence, perturbations in both hosts, virus and all host compartments interact. Different compartment residence times offer a possible interpretation of that coupled feedback, but no separate causal claim about a delay mechanism is proved here. The certified conclusion is that stable resident systems and positive invasion rates can coexist with an unstable interior equilibrium and, locally, stable oscillations.

One-host Hopf and enrichment behaviour already appears in Blath and Tóbiás (2023). Their Proposition 2.2 gives basic equilibrium spectral restrictions; Remark 2.3 and Section 3 discuss Hopf loss of stability and periodic behaviour, with rigorous results in special regimes alongside numerical investigations. The present contribution preserves two stable resident endemic systems and strict mutual invasion in a six-dimensional competing-host model, then certifies a bifurcation varying only dormancy. Thus the occurrence of a Hopf bifurcation in a host–virus model is established prior context, rather than the new claim.

The family resolves the specified deterministic question negatively with its resident-stability and invasion hypotheses preserved. The local bifurcation result identifies attracting periodic coexistence obtained by changing only dormancy. The family does not prove a periodic attractor for every $C\geq1/2$. The Hopf theorem does not give a global periodic continuation, an explicit neighbourhood radius, or uniform persistence for all positive trajectories. Neither theorem alters the source's proved stochastic early-invasion conclusions or settles its separate one-host global-convergence conjectures.

The supplement contains package-free JavaScript exact checks, symbolic rational verification, a standard-library integer-interval implementation, machine-readable certificates, and a runner recording commands, exit codes and source hashes. The accompanying provenance record specifies the source version and the fresh executions supporting this manuscript.

The proof is computer assisted, not proof-assistant formalisation. Its computational components are exposed for inspection. A separately transcribed exact witness check further guards against copying an erroneous Jacobian. Reproducing the proofs does not require the illustrative spectral plot.

A bounded literature search cannot prove historical priority. The contribution is the explicit equation-specific family and q-only Hopf result, using established Hurwitz, stable-manifold and normal-form methods.

## Appendix A. Exact resident characteristic polynomials

At $E_1$, the monic cubic has coefficients $1,c_1,c_2,c_3$:

$$
\begin{aligned}
c_1&=\frac{294C^2+17633C+895}{210C+150},\\
c_2&=\frac{353241C^2-105875C-1750}{1050C+750},\\
c_3&=\frac{11564}{25}C+7.
\end{aligned}
$$

The determinants are $c_1$, $c_1c_2-c_3$ and $c_3(c_1c_2-c_3)$. After $C=(t+1)/2$, their reduced numerators have strictly positive coefficients, with degrees $2,4,5$. Their denominators also have strictly positive coefficients.

At $E_2$, the monic quartic has coefficients $1,e_1,e_2,e_3,e_4$:

$$
\begin{aligned}
e_1&=\frac{4351200000C^2+362526500C+463091}{62160000C+50100},\\
e_2&=\frac{85544200000C^2+1109363500C+2375241}{207200000C+167000},\\
e_3&=\frac{231142800000C^2+1954099000C+5041897}{621600000C+501000},\\
e_4&=\frac{637}{4}C+\frac{1169}{2000}.
\end{aligned}
$$

Its determinants are $e_1$, $e_1e_2-e_3$, $e_3(e_1e_2-e_3)-e_1^2e_4$ and $e_4$ times the third determinant. The corresponding numerator degrees after the same substitution are $2,4,6,7$, again with strictly positive coefficients and positive-coefficient denominators. The family certificate prints these finite integer lists. They are computed from submatrices of the differentiated field.

At $C=1/2$, the exact resident viral densities are $V_1=1702/17$ and $V_2=38195/34589$. The latter is approximately $1.10425$ and differs from the resident viral density near the Hopf point, where $q$ has changed.

## Appendix B. The q-only characteristic coefficients

Set $L=5q+49$. Then

$$
\begin{aligned}
\alpha_1&=\frac{N_1}{300qL},&
\alpha_2&=\frac{N_2}{9000q^2L},&
\alpha_3&=\frac{N_3}{90000q^2L},\\
\alpha_4&=\frac{N_4}{2700000q^2L^2},&
\alpha_5&=\frac{-49N_5}{450000q^2L^2},&
\alpha_6&=\frac{833N_6}{600000q^2L},
\end{aligned}
$$

where

$$
\begin{aligned}
N_1={}&184865q^2+140287q+48020,\\
N_2={}&47850335q^3+32467841q^2+7087752q+1176490,\\
N_3={}&1429336425q^3+31184629q^2+1273247899q-39883011,\\
N_4={}&607919428125q^4-22001310797200q^3\\
&+29089900503650q^2-3765832488152q+269383385929,\\
N_5={}&14410621250q^4-19568133750q^3-25203074295q^2\\
&+21642947739q-1292021318,\\
N_6={}&(17020q-833)(17125q^2-16824q+833).
\end{aligned}
$$

The symbolic verifier proves these as rational-function identities, rather than from a finite sample alone. The integer polynomial script reconstructs $P$ directly: multiplying the Hurwitz matrix by $W=5400000q^2(5q+49)^2$ gives a degree-20 determinant numerator. Cancelling one factor $q$, two factors $5q+49$, and the integer content produces the primitive degree-17 polynomial. Exact endpoint and derivative evaluations then isolate its root.

## Appendix C. Bilinear form and executable proof

For $x,y\in\mathbb C^6$, define

$$
\begin{aligned}
S_x&=\sum_{k=1}^5x_k,\qquad S_y=\sum_{k=1}^5y_k,\\
A_{xy}&=x_1y_6+x_6y_1,\qquad B_{xy}=x_3y_6+x_6y_3.
\end{aligned}
$$

With $C=1/2,D=1$, the Hessian bilinear form is

$$
\begin{aligned}
B_1(x,y)&=-C(x_1S_y+y_1S_x)-A_{xy},\\
B_2(x,y)&=A_{xy},\\
B_3(x,y)&=-C(x_3S_y+y_3S_x)-B_{xy},\\
B_4(x,y)&=qB_{xy},\\
B_5(x,y)&=(1-q)B_{xy},\\
B_6(x,y)&=-A_{xy}-(1-q)B_{xy}.
\end{aligned}
$$

The symbolic verifier differentiates (1.1) to check this formula and that all third derivatives vanish.

\newpage

The portable runner fails on any unmet assertion in the following components.

| Component | Verified content |
| --- | --- |
| Rational witness | Equilibrium residuals, Jacobian, characteristic coefficients, trace, resident Hurwitz tests and full invasion indices at $C=1/2$ |
| Integer family determinants | All six determinant coefficient-sign certificates for $C\geq1/2$ |
| Exact q-root polynomial | Primitive $P$, endpoint signs and derivative exclusion |
| Rational q sample audit | Branch residuals and all characteristic coefficients at thirteen rational points |
| Symbolic identities | Full all-$C$ resident proof; q identities; imaginary factor modulo $P$; exact Hessian |
| Integer interval normal form | Stable quartic, frequency, eigenvectors, total crossing, $l_1$, and resident/invasion inequalities at the root |

The sample audit and spectral plot are supplementary. Every theorem-bearing claim has an exact or outward-interval certificate or an explicit analytical argument above.

## References

Blath, J., & Tóbiás, A. (2026). Emergence of microbial host dormancy during a persistent virus epidemic. *Journal of Mathematical Biology, 92*, Article 60. https://doi.org/10.1007/s00285-026-02372-8

Blath, J., & Tóbiás, A. (2023). Microbial virus epidemics in the presence of contact-mediated host dormancy. *ESAIM: Probability and Statistics, 27*, 174–220. https://doi.org/10.1051/ps/2022022

Gulbudak, H., & Weitz, J. S. (2016). A touch of sleep: Biophysical model of contact-mediated dormancy of archaea by viruses. *Proceedings of the Royal Society B: Biological Sciences, 283*(1839), 20161037. https://doi.org/10.1098/rspb.2016.1037

Kuznetsov, Y. A. (1998, November). *Explicit normal form coefficients for all codim 2 singularities of equilibria* [Lecture notes]. https://webspace.science.uu.nl/~kouzn101/INLDS/CODIM2.pdf

