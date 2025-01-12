Mass-balance model tailored to the Kaskawulsh Glacier, Yukon, Canada.

The distributed mass-balance model in this repository is adapted from Young et al. (2021). Several changes to the model workflow were introduced in Robinson et al. (2024), such as an annually updated surface elevation and a re-calibration of the model using a tuning scheme that includes distributed snowline observations. Major changes to the model include an improved treatment of supra-glacier debris based on site-specific measurements of ice ablation in the presence of debris and an improved elevation-dependent accumulation bias correction. See K. Robinson's MSc thesis for detailed descriptions of each model component. 

Robinson, K. (2024). Reconstructing a multi-decadal runoff record for a highly-glacierized catchment in Yukon, Canada. MSc thesis. https://summit.sfu.ca/item/38185

The model calculates the distributed surface mass balance $\dot{b}(x,y)$ as
\begin{equation}
    \label{massbalance}
    \dot{b}(x,y) = \dot{c}(x,y) - \dot{a}(x,y),
\end{equation}
where $\dot{c}(x,y)$ is the distributed surface accumulation and $\dot{a}(x,y)$ is the distributed surface ablation. Ablation is approximated as the difference between melt ($M$) and meltwater that refreezes to form superimposed ice ($R$). Melt is calculated using the enhanced temperature-index model of \citet{hock_distributed_1999}, 
\begin{equation}
\label{eq:TIM}
  M =
  \begin{cases}
    (MF + a_{\rm snow/ice}I)T & \text{if $T>0\degree\,C$} \\
    0 & \text{if $T\leqslant0 \degree\,C$},
  \end{cases}
\end{equation}
where T is air temperature and I is the potential direct clear-sky solar radiation (W\,m$^{-2}$). MF (m\,w.e.\,3hr$^{-1}$\,$\degree$C$^{-1}$), $a_{\rm snow}$ and $a_{\rm ice}$ (m\,w.e.\,3hr$^{-1}$\,$\degree$C$^{-1}$\,m$^2$\,W$^{-1}$) are, respectively, the melt factor and radiation factors for snow and ice that are empirically determined during the tuning process. 

The refreezing process is accounted for using a thermodynamic parameterization to estimate the total amount of liquid water (from snowmelt or rainfall) that can be retained by percolation and refreezing in the snowpack, hereafter referred to as the total potential retention mass P$_{\tau}$ \citep{janssens2000treatment}. This is equivalent to the maximum amount of superimposed ice that can form in a given year \citep{huybrechts1999dynamic}. For every hydrologic year (1 October--30 September)  in the study period, P$_{\tau}$ is approximated as a proportion ($P_r$) of the total annual precipitation in a given hydrological year ($P_{\rm annual}$): 
\begin{equation}
    \label{eq:refreeze}
    P_r = \frac{c}{L}\rm |min(T_{mean},0)|\frac{d}{P_{mean}},
\end{equation}
where $c$ is the specific heat capacity of ice, $L$ is the latent heat of fusion, $T_{\rm mean}$ is the local mean annual air temperature for a given hydrological year, $P_{\rm mean}$ is the mean annual precipitation over the whole study period (1980--2022) measured in m\,w.e., and $d$ is a prescribed thickness of the thermal active layer, set to 2\,m \citep{janssens2000treatment,young_imbalancing_2021}. The maximum allowable value of the retention fraction $P_r$ is 1, therefore the maximum possible  potential retention mass P$_{\tau}$ is equal to the annual precipitation ($P_{\rm annual}$), since
\begin{equation}
    \label{eq:Ptau}
    P_{\tau} = P_r\,P_{\rm annual}.
\end{equation}
While P$_{\tau}\,>\,0$, any melt that occurs is assumed to refreeze, therefore the maximum amount of refreezing that can occur is capped at P$_{\tau}$. Once the upper limit of P$_{\tau}$ has been reached, any additional snowmelt or rainfall is assumed to run off \citep{huybrechts1999dynamic,janssens2000treatment} until P$_{\tau}$ is renewed at the beginning of the next hydrological year. Therefore $R$ is related to the snowmelt ($M_{\rm snow}$) and P$_{\tau}$ in each gridcell and at each timestep by
\begin{equation}
\label{eq:refreezing}
  R =
  \begin{cases}
    M_{\rm snow} & \text{if $P_{\tau}$\,$\geq$\,$M_{\rm snow}$} \\
    P_{\tau} & \text{if 0\,$<$\,$P_{\tau}$\,$<$\,$M_{\rm snow}$} \\
    0 & \text{if $P_{\tau}\,=\,0$}.
  \end{cases}
\end{equation}
If the total refreezing $R\,<\,P_{\tau}$ for a given year, the remaining potential retention mass given by
\begin{equation}
    \label{eq:Ptau_remaining}
    P_{\tau\,\,\rm remaining} = P_{\tau} - R
\end{equation}
is carried over to the following hydrological year. 
