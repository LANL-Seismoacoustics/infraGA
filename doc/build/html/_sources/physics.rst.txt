.. _physics:

=====================================
Propagation Physics
=====================================

-------------------
Geometric Acoustics
-------------------

The propagation of acoustic energy can be described by a linear perturbation of the fluid mechanics equations.  The linear order continuity, Euler, and state equations for an inhomogeneous, moving medium have the forms,

    .. math:: 
         \frac{D \rho}{D t} + \rho \vec{\nabla} \cdot \vec{v}_0 + \vec{\nabla} \cdot \left( \rho_0 \vec{v} \right) & = 0,

         \frac{D \vec{v}}{Dt}  + \left( \vec{v} \cdot \vec{\nabla} \right) \vec{v}_0 & = -\frac{1}{\rho_0} \vec{\nabla} p + \rho \vec{\nabla} \frac{p_0}{\rho_0^2},

         \vec{v} \cdot \vec{\nabla} p_0 + \frac{D p}{D t} & = c^2 \left[ \vec{v} \cdot \vec{\nabla} \rho_0 + \frac{D \rho}{D t} \right] + \left( c^2 \right)^\prime \vec{v}_0 \cdot \vec{\nabla} \rho_0,


where subscript 0's denote ambient quantities that vary in space and those without subscripts denote acoustic perturbations that vary in space and time.  The approximation of geometric acoustics is constructed by expanding each variable with a spatially varying phase, :math:`e^{i k_0 \lambda \left( \vec{x} \right)}`, and Debye series, :math:`\sum{ \frac{\mathcal{P}_j (\vec{x})}{(i k_0)^j}}`.  The phase function, :math:`\psi \left( \vec{x} \right)`, is termed the Eikonal and its solution provides information about the deformation of surfaces of constant phase.  Expanding each linear variable,


    .. math::

        \begin{pmatrix} 
        p \\ \vec{v} \\ \rho \\ \left( c^2 \right)^\prime
        \end{pmatrix} =
        e^{i k_0 \lambda \left( \vec{x} \right)} \sum_{j = 0}^\infty{\frac{1}{\left( i k_0 \right)^j}
        \begin{pmatrix}
        \mathcal{P}_j \left( \vec{x} \right) \\
        \vec{\mathcal{V}}_j \left( \vec{x} \right) \\
        \mathcal{D}_j \left( \vec{x} \right) \\
        \mathcal{C}_j \left( \vec{x} \right)
        \end{pmatrix}}.

The linearized fluid mechanics equations can then be expressed as,

    .. math::

      & \sum_{j = 0}^\infty{\frac{1}{\left( i k_0 \right)^j} \Big\{ -i k_0 \mathcal{D}_j \left( c_0  - \vec{v}_0 \cdot \vec{\psi} \right) + \vec{v}_0 \cdot \vec{\nabla} \mathcal{D}_j + \mathcal{D}_j \vec{\nabla} \cdot \vec{v}_0}
      
      & \hspace{2.5in} + \rho_0 \vec{\nabla} \cdot \vec{\mathcal{V}}_j + \rho_0 i k_0 \vec{\mathcal{V}}_j \cdot \vec{\psi} + \vec{\mathcal{V}}_j \cdot \vec{\nabla} \rho_0 \Big\} = 0,

      &\sum_{j = 0}^\infty{\frac{1}{\left( i k_0 \right)^j} \left\{ -i k_0 \vec{\mathcal{V}}_j \left( c_0 - \vec{v}_0 \cdot \vec{\psi} \right) + \vec{v}_0 \cdot \vec{\nabla} \vec{\mathcal{V}}_j + \vec{\mathcal{V}}_j \cdot \vec{\nabla} \vec{v}_0 \right\}} 
      
      & \hspace{2.5in} = \sum_{j = 0}^\infty{\frac{1}{\left( i k_0 \right)^j}  \left\{ -  \frac{i k_0}{\rho_0}  \mathcal{P}_j   \vec{\psi} - \frac{1}{\rho_0} \vec{\nabla} \mathcal{P}_j  + \mathcal{D}_j \vec{\nabla} \frac{p_0}{\rho_0{}^2} \right\} },

      & \sum_{j = 0}^\infty{\frac{1}{\left( i k_0 \right)^j} \left\{ \vec{\mathcal{V}}_j \cdot \vec{\nabla} p_0 - i k_0 \mathcal{P}_j \left( c_0 -  \vec{v}_0 \cdot \vec{\psi} \right) + \vec{v}_0 \cdot \vec{\nabla} \mathcal{P}_j \right\} }
    
     & \hspace{1.0in} = \sum_{j = 0}^\infty{\frac{1}{\left( i k_0 \right)^j} \Big\{ c^2 \left[ \vec{\mathcal{V}}_j \cdot \vec{\nabla} \rho_0  -i k_0 \mathcal{D}_j \left( c_0 - \vec{v}_0 \cdot \vec{\psi} \right) + \vec{v}_0 \cdot \vec{\nabla} \mathcal{D}_j \right]} + \mathcal{C}_j \vec{v}_0 \cdot \vec{\nabla} \rho_0 \Big\},

where we've defined :math:`\vec{\psi} = \vec{\nabla} \lambda`.  Collecting terms in powers of :math:`k_0`, the leading order contributions require,

    .. math::
        \left( 1 - \frac{\vec{v}_0 \cdot \vec{\psi}}{c_0} \right) \mathcal{D}_0 & = \frac{\rho_0}{c_0} \vec{\mathcal{V}}_0 \cdot \vec{\psi},

        \left( 1 - \frac{\vec{v}_0 \cdot \vec{\psi}}{c_0} \right) \vec{\mathcal{V}}_0 & = \frac{1}{\rho_0 c_0} \mathcal{P}_0 \vec{\psi},
        
        \mathcal{P}_0 & = c^2 \mathcal{D}_0,

which can be combined to obtain the Eikonal Equation,

    .. math::
    
        \psi^2 = \frac{c_0^2}{c^2 \left( \vec{x} \right) } \left[ 1 - \frac{\vec{v}_0 \left( \vec{x} \right) \cdot \vec{\psi}}{c_0} \right]^2.


In addition to the ray path geometry defined by the Eikonal Equation, higher order terms in the expansion provide a means to estimate ray spreading and geometric attenuation.  Taking the terms in the expansion proportional to :math:`k_0`,

    .. math::
        
        & \left( 1 - \frac{\vec{v}_0 \cdot \vec{\psi}}{c_0} \right) \vec{\mathcal{V}}_1 - \frac{1}{\rho_0 c_0} \mathcal{P}_1 \vec{\psi}

        & \hspace{0.5in} =  \frac{1}{c_0} \left[ \vec{v}_0 \cdot\vec{\nabla}\vec{\mathcal{V}}_0 + \vec{\mathcal{V}}_0 \cdot\vec{\nabla}\vec{v}_0 + \frac{1}{\rho_0} \vec{\nabla} \mathcal{P}_0 - \frac{\mathcal{D}_0}{\rho_0{}^2}\vec{\nabla}p_0 \right] = \vec{b},

        &\left( 1 - \frac{\vec{v}_0 \cdot \vec{\psi}}{c_0} \right) \mathcal{D}_1 - \frac{\rho_0}{c_0} \vec{\psi} \cdot \vec{\mathcal{V}}_1 = \frac{1}{c_0}\vec{\nabla}\cdot \left( \mathcal{D}_0 \vec{v}_0 + \rho_0 \vec{\mathcal{V}}_0 \right) = b_1,

        &\mathcal{P}_1 - c^2 \mathcal{D}_1 = \frac{1}{c \psi} \left[ \vec{\mathcal{V}}_0 \cdot\vec{\nabla} p_0 + \vec{v}_0 \cdot\vec{\nabla}\mathcal{P}_0 - c^2 \vec{v}_0 \cdot\vec{\nabla} \mathcal{D}_0 - \frac{\mathcal{P}_0}{c^2} \vec{v}_0 \cdot \vec{\nabla} c^2 - c^2 \vec{\mathcal{V}}_0 \cdot\vec{\nabla}\rho_0 \right] = b_2,

Using the Eikonal Equation condition, these equations can be combined in a manner which goes to zero,

    .. math::

        \left\{
        \begin{matrix}
        \psi \vec{\mathcal{V}}_1 - \frac{1}{\rho_0 c_0} \mathcal{P}_1 \vec{\psi} = \vec{b} \\ \\
        \psi \mathcal{D}_1 - \frac{\rho_0}{c_0} \vec{\psi} \cdot \vec{\mathcal{V}}_1 = b_1 \\ \\
        \mathcal{P}_1 - c^2 \mathcal{D}_1 = b_2
        \end{matrix}
        \right. \quad \quad \rightarrow \quad
        \frac{c_0 \rho_0}{\psi} \vec{\psi} \cdot \vec{b} + c_0 c b_1 + \psi b_2 = 0.


Solving this condition leads to the transport equation,

    .. math::
        
        \vec{\nabla} \cdot \left( \mathcal{P}_0^2 \vec{c}_g \right) = \mathcal{P}_0^2 \vec{c}_g \cdot \vec{\nabla} \ln \left( \rho_0 c^3 \psi \right),

and the resulting amplitude term is defined in terms of the Jacobian, :math:`D \left( s, \vartheta, \varphi \right)` where :math`s`, :math:`\vartheta`, and :math:`\varphi` are the ray length, initial inclination angle, and initial azimuthal angle of the ray path respectively, that describes the coordinate transformation between Cartesian and ray coordinates, 


    .. math::

        \mathcal{P}_0 \left( s, \vartheta, \varphi \right) = \frac{1}{4 \pi} \sqrt{ \frac{\rho_0 \left( s \right) \psi \left( s \right) c^3 \left( s \right)}{\rho_0 \left( 0 \right) \psi \left( 0 \right) c^3 \left( 0 \right)} \frac{c_g \left( 0 \right) \,  \cos \vartheta}{c_g \left( s \right) D \left( s, \vartheta, \varphi \right) }}.


Spherical spreading at the source has been assumed so that :math:`\mathcal{P}_0 \left( s, \vartheta, \varphi \right)_{s \downarrow 0} = \frac{1}{4\pi s^2}` and :math:`D \left( s, \vartheta, \varphi \right)_{s \downarrow 0} = s^2 \cos \vartheta`, 


The Eikonal Equation can be used to define a Hamiltonian, :math:`  H \left( \vec{x}, \vec{\psi} \right) = 0 ` and the Hamilton-Jacobi relations used to define equations governing ray paths,

    .. math::
        
        \frac{\partial \vec{x}}{\partial \tau} & = \frac{\partial H}{\partial \vec{\psi}}, \quad \quad
        
        \frac{\partial \vec{\psi}}{\partial \tau} & = - \frac{\partial H}{\partial \vec{x}}.


For the methods in infraGA, the parameter :math:`\tau` is replaced by ray length, :math:`s`, such that :math:`\left\| d \vec{x} \right\| = ds`.  The transport coefficient depends on the Jacobian determinant which is defined by the variation between coordinate systems.  Denoting the initial launch inclination and azimuth as :math:`\vartheta` and :math:`\varphi`, respectively, one has,

    .. math::
        D \left( x, y, z; s, \vartheta, \varphi \right) = \left\| \frac{\partial \left( x, y, z \right)}{\partial \left( s, \vartheta, \varphi \right)} \right\| = \textbf{det} \begin{pmatrix}
        \frac{\partial x}{\partial s} && \frac{\partial x}{\partial \vartheta} && \frac{\partial x}{\partial \varphi} \\ \\
        \frac{\partial y}{\partial s} && \frac{\partial y}{\partial \vartheta} && \frac{\partial y}{\partial \varphi} \\ \\
        \frac{\partial z}{\partial s} && \frac{\partial z}{\partial \vartheta} && \frac{\partial z}{\partial \varphi}
        \end{pmatrix}

As detailed in Blom \& Waxler (2012), the :math:`s` derivatives can be defined directly from the Eikonal Equation condition, but the :math:`\vartheta` and :math:`\varphi` derivatives require the introduction of auxiliary parameters, :math:`\mathcal{X}^{(\vartheta)} = \frac{\partial x}{\partial \vartheta}`, :math:`\mathcal{X}^{(\varphi)} = \frac{\partial x}{\partial \varphi}` with similar parameters defined for :math:`\mathcal{Y}` and :math:`\mathcal{Z}`.  The governing equations for the auxiliary parameters are defined by taking launch angle derivatives of the governing spatial and Eikonal differential equations, for example,

    .. math::
        \frac{\partial \mathcal{X}^{(\vartheta)}}{\partial s} & = \frac{\partial}{\partial \vartheta} \frac{\partial x}{\partial s}, \\ 
        \frac{\partial \Phi_x^{(\vartheta)}}{\partial s} & = \frac{\partial}{\partial \vartheta} \frac{\partial \psi_x}{\partial s}

This increases the number of coupled equations needed for ray computation by a factor of 2 for 2-dimensional simulations and by a factor of 3 for 3-dimensional and spherical geometry simulations in which variations with respect to both :math:`\vartheta` and :math:`\varphi` must be considered.

--------------------------------------
Two- and Three-Dimensional Propagation
--------------------------------------

In the case of the effective sound speed approximation, one re-defines :math:`c \rightarrow c + \vec{v}_0 \cdot \hat{\psi}_\perp` (adding the wind in the direction of propagation to the adiabatic sound speed) and :math:`\vec{v} = 0` in the relations.  This reduces the Eikonal to, :math:`\psi^2 = \frac{c_0^2}{c^2}`, and the propagation relations become simply,

    .. math::

        \frac{\partial \vec{x}}{\partial s} = \frac{c_0}{c} \vec{\psi} , \quad \quad \frac{\partial \psi_j}{\partial s} = - \frac{c_0}{c^2} \frac{\partial c}{\partial x_j} ,

For three-dimensional propagation simulations, the differential equations describing the geometric ray paths in an arbitrary moving medium can be found from the unmodified Eikonal derived above,

    .. math::

        \frac{\partial \vec{x}}{\partial s} & = \frac{\vec{c}_g}{c_g}, \quad \vec{c}_g = c \frac{\vec{\psi}}{\psi} + \vec{v}_0
        
        \frac{\partial \psi_j}{\partial s} & = - \frac{1}{c_g} \left[ \psi \frac{\partial c}{\partial x_j} + \vec{\psi} \cdot \frac{\partial \vec{v}_0}{\partial x_j} \right].

where :math:`\vec{c}_g` is the group velocity of energy along the ray paths.  See Blom & Waxler 2012 and 2017 for full discussion of the Cartesian ray tracing development and inclusion of auxiliary parameters.

------------------------------
Spherical Geometry Propagation
------------------------------

The eikonal solution in spherical coordinates requires geometric corrections to the scaling of :math:`d\vec{x}` as well as additional terms in the :math:`\frac{\partial \psi_j}{\partial s}` relations to preserve the Eikonal vector direction as unit vectors vary in space,

    .. math::
        \frac{\partial u_j}{\partial s} & = \mathcal{G}_j\frac{c_{g,j}}{c_g} , \quad \quad c_{g,j} = c \frac{\psi_j}{\psi} + v_j,

        \frac{\partial \psi_j}{\partial s} & = - \frac{\mathcal{G}_j}{c_g} \left( \psi \frac{\partial c}{\partial u_j} + \sum_k{ \psi_k \frac{\partial v_k}{\partial u_j}} + \mathcal{T}_j \right),

where geometric scaling coefficients and corrective terms for the spatial variability of the spherical coordinate unit vectors produce,

    .. math::
        \mathcal{G}_r & = 1,  \quad 
        \mathcal{G}_\theta = \frac{1}{r}, \quad
        \mathcal{G}_\phi = \frac{1}{r \cos \theta},
        
        \mathcal{T}_r & =  \frac{1}{r} \left( \psi_\theta c_{g,\theta} + \psi_\phi  c_{g,\phi} \right),

        \mathcal{T}_\theta & =  \left( \psi_r v_\theta - \psi_\theta v_r \right) - \left( \psi_r c_{g,\theta} - \psi_\phi  c_{g,\phi}  \tan \theta \right),

        \mathcal{T}_\phi & =  \left( \psi_r v_\phi - \psi_\phi v_r \right) \cos \theta

        & \hspace{0.5in} + \left( \psi_\theta v_\phi - \psi_\phi v_\theta \right) \sin \theta

        & \hspace{1.0in} - c_{g,\phi} \left( \psi_r \cos \theta + \psi_\theta \sin \theta \right).


This is the equation set used in the *infraga-sph* methods and the full derivation of these relations is included in Blom, 2019.

----------------
Eigenray Methods
----------------


Eigenrays are identified using the auxiliary parameters defined in order to compute the Jacobian components needed to calculate geometric spreading.  Considering the arrival location of a ray path in 3D,


    .. math::
        
        x_0 \left( \vartheta + \delta \vartheta, \varphi + \delta \varphi \right) = x_0 \left( \vartheta, \varphi \right) + \frac{\partial x_0}{\partial \vartheta} \delta \vartheta + \frac{\partial x_0}{\partial \varphi} \delta \varphi + O \left( \delta^2 \right),  \\
        
        y_0 \left( \vartheta + \delta \vartheta, \varphi + \delta \varphi \right) = y_0 \left( \vartheta, \varphi \right) + \frac{\partial y_0}{\partial \vartheta} \delta \vartheta + \frac{\partial y_0}{\partial \varphi} \delta \varphi + O \left( \delta^2 \right),

The arrival location shifts with respect to launch angle variations (:math:`\frac{\partial x_0}{\partial \vartheta}`, :math:`\frac{\partial x_0}{\partial \varphi}`, etc.) can be defined from the auxiliary parameters introduced to solve the Transport equation.  This can be written more compactly as,

    .. math::
        \begin{pmatrix}
        \delta x_0 \\
        \delta y_0  
        \end{pmatrix} = 
        \begin{pmatrix}
        \frac{\partial x_0}{\partial \vartheta} & \frac{\partial x_0}{\partial \varphi} \\
        \frac{\partial y_0}{\partial \vartheta} & \frac{\partial y_0}{\partial \varphi}
        \end{pmatrix}
        \begin{pmatrix}
        \delta \vartheta \\ \delta \varphi
        \end{pmatrix},

or simply,

    .. math::

        \delta \vec{x}_0 = \boldsymbol{\mathcal{D}}_0 \, \delta \vec{\vartheta}.

From this linear approximation, a Levenberg-Marquardt algorithm can be constructed,

    .. math::
        
        \delta \vec{\vartheta} = \left( \boldsymbol{\mathcal{D}}_0 + \lambda  \text{ diag} \left(\boldsymbol{\mathcal{D}}_0 \right) \right)^{-1} \delta \vec{x}_0,


that will identify the changes in ray launch angles, :math:`\delta \vec{\vartheta}` needed to shift the arrival location by some distance, :math:`\delta \vec{x}_0`.  This algorithm is utilized as a stand alone method in *-eig_direct* and as the precision search step in *-eig_search* where a preliminary inclination/range search is used to identify initial solutions near eigenrays.  See Blom & Waxler 2017 for a full discussion of the eigenray methods.

--------------------
Terrain Interactions
--------------------

The reflection conditions for including topography are modified so that the eikonal vector components along the ground surface are conserved and the normal component to the ground changes sign,


    .. math::

        \vec{\psi}_\text{refl} \cdot \hat{n}_\text{grnd} = - \vec{\psi}_0 \cdot \hat{n}_\text{grnd}, \quad
        \vec{\psi}_\text{refl} \times \hat{n}_\text{grnd} = \vec{\psi}_0 \times \hat{n}_\text{grnd}.
    
where :math:`\vec{\psi}_0` denotes the incident eikonal vector.  The resulting reflection conditions can then be defined by relating the ground normal to the derivative of the functional ground surface specification, :math:`z_g \left( x, y \right)`,


    .. math::

        \psi_x \left( s_0 + 0^+, \vartheta, \varphi \right) & = \mathcal{C}_1^{(x)}  \psi_{x,0} + \mathcal{C}_2^{(x)} \left( \psi_{z,0} - \psi_{y,0} \frac{\partial z_g}{\partial y} \right), \\
        \psi_y \left( s_0 + 0^+, \vartheta, \varphi \right) & = \mathcal{C}_1^{(y)}  \psi_{y,0} + \mathcal{C}_2^{(y)} \left( \psi_{z,0} - \psi_{x,0} \frac{\partial z_g}{\partial x} \right), \\        
        \psi_z \left( s_0 + 0^+, \vartheta, \varphi \right) & = -\mathcal{C}_1^{(z)} \psi_{z, 0} + \mathcal{C}_2^{(x)} \psi_{x,0} + \mathcal{C}_2^{(y)} \psi_{y,0}.

where,

    .. math::
        
        \mathcal{C}_1^{(x)} & = \frac{1 - \left(\frac{\partial z_g}{\partial x} \right)^2 + \left(\frac{\partial z_g}{\partial y} \right)^2}{1 + \left(\frac{\partial z_g}{\partial x} \right)^2 + \left(\frac{\partial z_g}{\partial y} \right)^2}, & \quad \quad
        \mathcal{C}_2^{(x)} &= \frac{2 \frac{\partial z_g}{\partial x}}{1 + \left(\frac{\partial z_g}{\partial x} \right)^2 +\left(\frac{\partial z_g}{\partial y} \right)^2}, \\
        \mathcal{C}_1^{(y)} & =\frac{1 + \left(\frac{\partial z_g}{\partial x} \right)^2 - \left(\frac{\partial z_g}{\partial y} \right)^2}{1 + \left(\frac{\partial z_g}{\partial x} \right)^2 + \left(\frac{\partial z_g}{\partial y} \right)^2}, & \quad \quad 
        \mathcal{C}_2^{(y)} & = \frac{2 \frac{\partial z_g}{\partial y}}{1 + \left(\frac{\partial z_g}{\partial x} \right)^2 +\left(\frac{\partial z_g}{\partial y} \right)^2}, \\
        \mathcal{C}_1^{(z)} & = \frac{1 - \left(\frac{\partial z_g}{\partial x} \right)^2 - \left(\frac{\partial z_g}{\partial y} \right)^2}{1 + \left(\frac{\partial z_g}{\partial x} \right)^2 + \left(\frac{\partial z_g}{\partial y} \right)^2}.


See Blom 2020 for a full derivation and discussion of the reflection conditions.


------------------------------------
Weakly Non-Linear Waveform Evolution
------------------------------------

The waveform evolution is computed using the methods developed by Lonzaga \textit{et al.} (2015) using a Heun's solver (RK2).  The equations being solved are,


    .. math::
        \frac{\partial u}{\partial s} = \tilde{\beta} u \frac{\partial u}{\partial \tau}, \quad 
        \tilde{\beta} \left( s\right) & = \beta \frac{p_0}{\rho_0 c_0^2} \frac{\psi_0 c_0}{c_{g,0} c_\text{src}} \sqrt{ \frac{D_0 \rho_0 c_{g,0}^3}{D \rho c_g^3} \frac{c \psi^3}{c_0 \psi_0^3}}, \\
        u \left( s, \tau \right) & = \frac{p \left( s, \tau \right)}{p_\text{ref}} \sqrt{ \frac{\rho_0 c_0^3 \psi_0}{\rho c^3 \psi} \frac{c_{g} D}{c_{g,0} D_0}},

where subscript zeros denote evaluation at some reference point, :math:`s = s_0`, along the ray path.  The variable step size in the solver is defined as :math:`ds = ds_0 /  \left( \pi \tilde{\beta} \left( s \right) \text{max} \left(\mathcal{U} \left( s, f \right) \right) \right)`, where :math:`\mathcal{U} \left( s, f \right)` is the FFT of :math:`u \left(s, t \right)` along the ray path and :math:`ds_0` is defined in the code as \verb=wvfrm_ds=.  

In cases for which little energy is expected hear the Nyquist frequency, a value of :math:`ds_0 \sim1.0` can be used; however, for source waveforms with high frequency content (e.g., a blast wave) or propagation paths extending into the upper atmosphere where rarefaction leads to strong relatively strong non-linear effects and generation of high frequency energy, a value of :math:`\sim0.1` might be required.  Best practice is to vary the value of \verb=wvfrm_ds= to be sure your analysis has converged.  See Lonzaga et al., 2015 and Blom & Waxler, 2021 for a full description of the Burgers equation methods.

The blastwave/impulse source available in the software has the form,

    .. math::
         p \left( t; p_0, t_0, \alpha \right) = \left\{
        \begin{matrix}
        \frac{p_0}{\mathcal{C}} x^\alpha \left( 1 - \frac{x}{1 + \alpha} \right) e^{-x} &  x \geq 0 \\ \\
        0									& x < 0
        \end{matrix}
        \right., \quad \quad x = \frac{t}{t_0}.

for peak overpressure, :math:`p_0`, positive phase duration, :math:`t_0`, and shaping parameter, :math:`\alpha`.  This waveform was introduced by Waxler & Assink (2018) an improved source model for waveform simulations as it avoids the symmetry and narrow-banded limitations of a Gaussian enveloped sinusoid.  Interestingly, when the shaping parameter, :math:`\alpha`, approaches zero, the impulse becomes the Friedlander (1946) blastwave.

