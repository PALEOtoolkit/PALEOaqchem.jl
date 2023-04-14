var documenterSearchIndex = {"docs":
[{"location":"PALEOaqchem functions/#PALEOaqchem-functions","page":"PALEOaqchem functions","title":"PALEOaqchem functions","text":"","category":"section"},{"location":"PALEOaqchem functions/","page":"PALEOaqchem functions","title":"PALEOaqchem functions","text":"CurrentModule = PALEOaqchem","category":"page"},{"location":"PALEOaqchem functions/","page":"PALEOaqchem functions","title":"PALEOaqchem functions","text":"O2AlkUptakeRemin\n\nStoichPars","category":"page"},{"location":"PALEOaqchem functions/#PALEOaqchem.O2AlkUptakeRemin","page":"PALEOaqchem functions","title":"PALEOaqchem.O2AlkUptakeRemin","text":"O2AlkUptakeRemin(Corg, (NO3, TNH3, Ngas), TPO4, Ccarb; rO2Corg=1) -> (O2, Alk)\n\nOxygen and alkalinity assimilated for production (or released by remineralisation) of particulate matter with specified Corg, Ccarb from specified nitrate NO3, total ammonia TNH3, gaseous N, total phosphate TPO4\n\nNB sign: for Corg +ve, O2 is -ve (ie quantity to subtract from tracer sms for production, or add for remineralisation)\n\n\n\n\n\n","category":"function"},{"location":"PALEOaqchem functions/#PALEOaqchem.StoichPars","page":"PALEOaqchem functions","title":"PALEOaqchem.StoichPars","text":"StoichPars() -> PB.ParametersTuple\n\nParameters defining particulate organic matter stoichiometry\n\n\n\n\n\n","category":"function"},{"location":"References/#References","page":"References","title":"References","text":"","category":"section"},{"location":"References/","page":"References","title":"References","text":"","category":"page"},{"location":"indexpage/#Index","page":"Index","title":"Index","text":"","category":"section"},{"location":"indexpage/","page":"Index","title":"Index","text":"","category":"page"},{"location":"PALEOaqchem Reactions/#PALEOaqchem-Reactions","page":"PALEOaqchem Reactions","title":"PALEOaqchem Reactions","text":"","category":"section"},{"location":"PALEOaqchem Reactions/","page":"PALEOaqchem Reactions","title":"PALEOaqchem Reactions","text":"CurrentModule = PALEOaqchem","category":"page"},{"location":"PALEOaqchem Reactions/#Particulate-fluxes","page":"PALEOaqchem Reactions","title":"Particulate fluxes","text":"","category":"section"},{"location":"PALEOaqchem Reactions/","page":"PALEOaqchem Reactions","title":"PALEOaqchem Reactions","text":"Particle.ReactionParticleDecay\nParticle.ReactionFluxToComponents","category":"page"},{"location":"PALEOaqchem Reactions/#PALEOaqchem.Particle.ReactionParticleDecay","page":"PALEOaqchem Reactions","title":"PALEOaqchem.Particle.ReactionParticleDecay","text":"ReactionParticleDecay\n\nDecay (eg organic matter to remineralization) at rate 1.0/decay_timescale of eg an organic matter dissolved/particulate phase in Reservoir Particle, to decayflux. Particle may have isotope_type. The Reservoir Particle may have the :vsink attribute set to represent a marine sinking particulate phase.\n\nParameters\n\ndecay_timescale[Float64]=0.5  (yr), default_value=0.5, description=\"particle decay timescale\"\ndecay_threshold[Float64]=-Inf  (mol m-3), default_value=-Inf, description=\"particle decay concentration below which decay stops\"\nfield_data[DataType]=PALEOboxes.ScalarData, default_value=PALEOboxes.ScalarData, allowed_values=Type[PALEOboxes.ScalarData, PALEOboxes.IsotopeLinear], description=\"disable / enable isotopes and specify isotope type\"\n\nMethods and Variables\n\ndo_particle_decay\nParticle  (mol), VT_ReactDependency, description=\"Particle amount\"\nParticle_sms  (mol yr-1), VT_ReactContributor, description=\"Particle source-sink\"\ndecayflux  (mol yr-1), VT_ReactContributor, description=\"Particle decay flux\"\n\n\n\n\n\n","category":"type"},{"location":"PALEOaqchem Reactions/#PALEOaqchem.Particle.ReactionFluxToComponents","page":"PALEOaqchem Reactions","title":"PALEOaqchem.Particle.ReactionFluxToComponents","text":"ReactionFluxToComponents\n\nDistribute a single input_flux (eg an organic matter flux) to output_fluxnames components with stoichiometry output_fluxstoich. input_flux may have an isotope type (set by field_data) in which case the isotopic composition will be sent to (usually one) output_fluxname with ::Isotope suffix.\n\nParameters\n\noutputflux_prefix[String]=\"\", default_value=\"\", description=\"Prefix for output flux names\"\noutputflux_names[Vector{String}]=[\"Corg\", \"N\", \"P\"], default_value=[\"Corg\", \"N\", \"P\"], description=\"Suffixes for output flux names.  Use ::Isotope suffix to identify a flux with 'isotope_type'\"\noutputflux_stoich[Vector{Float64}]=[106.0, 16.0, 1.0], default_value=[106.0, 16.0, 1.0], description=\"Stoichiometry for output fluxes relative to input flux\"\nfield_data[DataType]=PALEOboxes.ScalarData, default_value=PALEOboxes.ScalarData, allowed_values=Type[PALEOboxes.ScalarData, PALEOboxes.IsotopeLinear], description=\"disable / enable isotopes and specify isotope type\"\n\nMethods and Variables for default Parameters\n\ndo_flux_to_components\ninputflux  (mol m-3), VT_ReactTarget, description=\"input flux\"\n[Corg]  (mol yr-1), VT_ReactContributor, description=\"flux Corg\"\n[N]  (mol yr-1), VT_ReactContributor, description=\"flux N\"\n[P]  (mol yr-1), VT_ReactContributor, description=\"flux P\"\n\n\n\n\n\n","category":"type"},{"location":"PALEOaqchem Reactions/#Remin","page":"PALEOaqchem Reactions","title":"Remin","text":"","category":"section"},{"location":"PALEOaqchem Reactions/","page":"PALEOaqchem Reactions","title":"PALEOaqchem Reactions","text":"Remin.ReactionReminPonly\nRemin.ReactionReminO2\nRemin.ReactionReminO2_SO4\nRemin.ReactionReminO2_SO4_CH4","category":"page"},{"location":"PALEOaqchem Reactions/#PALEOaqchem.Remin.ReactionReminPonly","page":"PALEOaqchem Reactions","title":"PALEOaqchem.Remin.ReactionReminPonly","text":"ReactionReminPonly\n\nOrganic particulate matter remineralization (no oxidant use)\n\nParameters\n\nCIsotope[external, DataType]=PALEOboxes.ScalarData, default_value=PALEOboxes.ScalarData, allowed_values=Type[PALEOboxes.ScalarData, PALEOboxes.IsotopeLinear], description=\"disable / enable carbon isotopes and specify isotope type\"\n\nMethods and Variables for default Parameters\n\ndo_remin_Ponly\nremin_P  (mol yr-1), VT_ReactTarget, description=\"flux P\"\nremin_N  (mol yr-1), VT_ReactTarget, description=\"flux N\"\nremin_Corg  (mol yr-1), VT_ReactTarget, description=\"flux Corg\"\nremin_Ccarb  (mol yr-1), VT_ReactTarget, description=\"flux Ccarb\"\nsoluteflux_P  (mol yr-1), VT_ReactContributor, description=\"flux P\"\n[soluteflux_TNH3]  (mol yr-1), VT_ReactContributor, description=\"flux TNH3\"\n[soluteflux_DIC]  (mol yr-1), VT_ReactContributor, description=\"flux DIC\"\n[soluteflux_TAlk]  (mol yr-1), VT_ReactContributor, description=\"flux TAlk\"\n\n\n\n\n\n","category":"type"},{"location":"PALEOaqchem Reactions/#PALEOaqchem.Remin.ReactionReminO2","page":"PALEOaqchem Reactions","title":"PALEOaqchem.Remin.ReactionReminO2","text":"ReactionReminO2\n\nOrganic particulate matter remineralization (O2 oxidant only)\n\nParameters\n\nCIsotope[external, DataType]=PALEOboxes.ScalarData, default_value=PALEOboxes.ScalarData, allowed_values=Type[PALEOboxes.ScalarData, PALEOboxes.IsotopeLinear], description=\"disable / enable carbon isotopes and specify isotope type\"\n\nMethods and Variables for default Parameters\n\ndo_remin_O2\nremin_P  (mol yr-1), VT_ReactTarget, description=\"flux P\"\nremin_N  (mol yr-1), VT_ReactTarget, description=\"flux N\"\nremin_Corg  (mol yr-1), VT_ReactTarget, description=\"flux Corg\"\nremin_Ccarb  (mol yr-1), VT_ReactTarget, description=\"flux Ccarb\"\nsoluteflux_P  (mol yr-1), VT_ReactContributor, description=\"flux P\"\n[soluteflux_TNH3]  (mol yr-1), VT_ReactContributor, description=\"flux TNH3\"\n[soluteflux_DIC]  (mol yr-1), VT_ReactContributor, description=\"flux DIC\"\n[soluteflux_TAlk]  (mol yr-1), VT_ReactContributor, description=\"flux TAlk\"\nreminOrgOxO2  (mol O2 yr-1), VT_ReactProperty, description=\"oxygen consumption (-ve) by remineralization\"\nRateStoich_reminOrgOxO2\nreminOrgOxO2  (mol O2 yr-1), VT_ReactDependency, description=\"oxygen consumption (-ve) by remineralization\"\n[soluteflux_O2]  (), VT_ReactContributor, description=\"generated by RateStoich rate=reminOrgOxO2\"\ntotals\nreminOrgOxO2  (mol O2 yr-1), VT_ReactDependency, description=\"oxygen consumption (-ve) by remineralization\"\nreminOrgOxO2_total  (mol O2 yr-1), VT_ReactProperty, description=\"total oxygen consumption (-ve) by remineralization\"\n\n\n\n\n\n","category":"type"},{"location":"PALEOaqchem Reactions/#PALEOaqchem.Remin.ReactionReminO2_SO4","page":"PALEOaqchem Reactions","title":"PALEOaqchem.Remin.ReactionReminO2_SO4","text":"ReactionReminO2_SO4\n\nOrganic particulate matter remineralization (O2, SO4 oxidants)\n\nParameters\n\noxreminlimit[Float64]=0.008  (mol m-3), default_value=0.008, description=\"oxygen concentration below which use of O2 for remineralisation is inhibited\"\nCIsotope[external, DataType]=PALEOboxes.ScalarData, default_value=PALEOboxes.ScalarData, allowed_values=Type[PALEOboxes.ScalarData, PALEOboxes.IsotopeLinear], description=\"disable / enable carbon isotopes and specify isotope type\"\nSIsotope[external, DataType]=PALEOboxes.ScalarData, default_value=PALEOboxes.ScalarData, allowed_values=Type[PALEOboxes.ScalarData, PALEOboxes.IsotopeLinear], description=\"disable / enable sulphur isotopes and specify isotope type\"\n\nMethods and Variables for default Parameters\n\ndo_remin_O2_SO4\nremin_P  (mol yr-1), VT_ReactTarget, description=\"flux P\"\nremin_N  (mol yr-1), VT_ReactTarget, description=\"flux N\"\nremin_Corg  (mol yr-1), VT_ReactTarget, description=\"flux Corg\"\nremin_Ccarb  (mol yr-1), VT_ReactTarget, description=\"flux Ccarb\"\nsoluteflux_P  (mol yr-1), VT_ReactContributor, description=\"flux P\"\n[soluteflux_TNH3]  (mol yr-1), VT_ReactContributor, description=\"flux TNH3\"\n[soluteflux_DIC]  (mol yr-1), VT_ReactContributor, description=\"flux DIC\"\n[soluteflux_TAlk]  (mol yr-1), VT_ReactContributor, description=\"flux TAlk\"\nO2_conc  (mol m-3), VT_ReactDependency, description=\"O2 concentration\"\nreminOrgOxO2  (mol O2 yr-1), VT_ReactProperty, description=\"oxygen consumption (-ve) by remineralization\"\nreminOrgOxSO4  (mol O2eq yr-1), VT_ReactProperty, description=\"2 * sulphate consumption (-ve) by remineralization\"\nRateStoich_reminOrgOxO2\nreminOrgOxO2  (mol O2 yr-1), VT_ReactDependency, description=\"oxygen consumption (-ve) by remineralization\"\n[soluteflux_O2]  (), VT_ReactContributor, description=\"generated by RateStoich rate=reminOrgOxO2\"\nRateStoich_reminOrgOxSO4\nreminOrgOxSO4  (mol O2eq yr-1), VT_ReactDependency, description=\"2 * sulphate consumption (-ve) by remineralization\"\n[soluteflux_SO4]  (), VT_ReactContributor, description=\"generated by RateStoich rate=reminOrgOxSO4\"\n[soluteflux_H2S]  (), VT_ReactContributor, description=\"generated by RateStoich rate=reminOrgOxSO4\"\n[soluteflux_TAlk]  (), VT_ReactContributor, description=\"generated by RateStoich rate=reminOrgOxSO4\"\ntotals\nreminOrgOxSO4  (mol O2eq yr-1), VT_ReactDependency, description=\"2 * sulphate consumption (-ve) by remineralization\"\nreminOrgOxO2  (mol O2 yr-1), VT_ReactDependency, description=\"oxygen consumption (-ve) by remineralization\"\nreminOrgOxSO4_total  (mol O2eq yr-1), VT_ReactProperty, description=\"total 2 * sulphate consumption (-ve) by remineralization\"\nreminOrgOxO2_total  (mol O2 yr-1), VT_ReactProperty, description=\"total oxygen consumption (-ve) by remineralization\"\n\n\n\n\n\n","category":"type"},{"location":"PALEOaqchem Reactions/#PALEOaqchem.Remin.ReactionReminO2_SO4_CH4","page":"PALEOaqchem Reactions","title":"PALEOaqchem.Remin.ReactionReminO2_SO4_CH4","text":"ReactionReminO2_SO4_CH4\n\nOrganic particulate matter remineralization (O2, SO4 oxidants, remaining Corg to CH4)\n\nParameters\n\noxreminlimit[Float64]=0.008  (mol m-3), default_value=0.008, description=\"oxygen concentration below which use of O2 for remineralisation is inhibited\"\nSO4reminlimit[Float64]=1.0  (mol m-3), default_value=1.0, description=\"sulphate concentration below which use of SO4 for remineralisation is inhibited\"\nCIsotope[external, DataType]=PALEOboxes.ScalarData, default_value=PALEOboxes.ScalarData, allowed_values=Type[PALEOboxes.ScalarData, PALEOboxes.IsotopeLinear], description=\"disable / enable carbon isotopes and specify isotope type\"\nSIsotope[external, DataType]=PALEOboxes.ScalarData, default_value=PALEOboxes.ScalarData, allowed_values=Type[PALEOboxes.ScalarData, PALEOboxes.IsotopeLinear], description=\"disable / enable sulphur isotopes and specify isotope type\"\n\nMethods and Variables for default Parameters\n\ndo_remin_O2_SO4_CH4\nremin_P  (mol yr-1), VT_ReactTarget, description=\"flux P\"\nremin_N  (mol yr-1), VT_ReactTarget, description=\"flux N\"\nremin_Corg  (mol yr-1), VT_ReactTarget, description=\"flux Corg\"\nremin_Ccarb  (mol yr-1), VT_ReactTarget, description=\"flux Ccarb\"\nsoluteflux_P  (mol yr-1), VT_ReactContributor, description=\"flux P\"\n[soluteflux_TNH3]  (mol yr-1), VT_ReactContributor, description=\"flux TNH3\"\n[soluteflux_DIC]  (mol yr-1), VT_ReactContributor, description=\"flux DIC\"\n[soluteflux_TAlk]  (mol yr-1), VT_ReactContributor, description=\"flux TAlk\"\nO2_conc  (mol m-3), VT_ReactDependency, description=\"O2 concentration\"\nSO4_conc  (mol m-3), VT_ReactDependency, description=\"SO4 concentration\"\nreminOrgOxO2  (mol O2 yr-1), VT_ReactProperty, description=\"oxygen consumption (-ve) by remineralization\"\nreminOrgOxSO4  (mol O2eq yr-1), VT_ReactProperty, description=\"2 * sulphate consumption (-ve) by remineralization\"\nreminOrgOxCH4  (mol O2eq yr-1), VT_ReactProperty, description=\"2 * DIC -> methane (-ve) by remineralization\"\nRateStoich_reminOrgOxO2\nreminOrgOxO2  (mol O2 yr-1), VT_ReactDependency, description=\"oxygen consumption (-ve) by remineralization\"\n[soluteflux_O2]  (), VT_ReactContributor, description=\"generated by RateStoich rate=reminOrgOxO2\"\nRateStoich_reminOrgOxSO4\nreminOrgOxSO4  (mol O2eq yr-1), VT_ReactDependency, description=\"2 * sulphate consumption (-ve) by remineralization\"\n[soluteflux_SO4]  (), VT_ReactContributor, description=\"generated by RateStoich rate=reminOrgOxSO4\"\n[soluteflux_H2S]  (), VT_ReactContributor, description=\"generated by RateStoich rate=reminOrgOxSO4\"\n[soluteflux_TAlk]  (), VT_ReactContributor, description=\"generated by RateStoich rate=reminOrgOxSO4\"\nRateStoich_reminOrgOxCH4\nreminOrgOxCH4  (mol O2eq yr-1), VT_ReactDependency, description=\"2 * DIC -> methane (-ve) by remineralization\"\n[soluteflux_DIC]  (), VT_ReactContributor, description=\"generated by RateStoich rate=reminOrgOxCH4\"\n[soluteflux_CH4]  (), VT_ReactContributor, description=\"generated by RateStoich rate=reminOrgOxCH4\"\ntotals\nreminOrgOxSO4  (mol O2eq yr-1), VT_ReactDependency, description=\"2 * sulphate consumption (-ve) by remineralization\"\nreminOrgOxO2  (mol O2 yr-1), VT_ReactDependency, description=\"oxygen consumption (-ve) by remineralization\"\nreminOrgOxCH4  (mol O2eq yr-1), VT_ReactDependency, description=\"2 * DIC -> methane (-ve) by remineralization\"\nreminOrgOxSO4_total  (mol O2eq yr-1), VT_ReactProperty, description=\"total 2 * sulphate consumption (-ve) by remineralization\"\nreminOrgOxO2_total  (mol O2 yr-1), VT_ReactProperty, description=\"total oxygen consumption (-ve) by remineralization\"\nreminOrgOxCH4_total  (mol O2eq yr-1), VT_ReactProperty, description=\"total 2 * DIC -> methane (-ve) by remineralization\"\n\n\n\n\n\n","category":"type"},{"location":"PALEOaqchem Reactions/#Secondary-redox","page":"PALEOaqchem Reactions","title":"Secondary redox","text":"","category":"section"},{"location":"PALEOaqchem Reactions/","page":"PALEOaqchem Reactions","title":"PALEOaqchem Reactions","text":"SecondaryRedox.ReactionRedoxH2S_O2\nSecondaryRedox.ReactionRedoxCH4_O2\nSecondaryRedox.ReactionRedoxCH4_SO4","category":"page"},{"location":"PALEOaqchem Reactions/#PALEOaqchem.SecondaryRedox.ReactionRedoxH2S_O2","page":"PALEOaqchem Reactions","title":"PALEOaqchem.SecondaryRedox.ReactionRedoxH2S_O2","text":"ReactionRedoxH2S_O2\n\nSulphide oxidation by oxygen.\n\nRate R_H2S_O2 Units Ref Notes\n1.6e5 (mol l-1)-1 yr-1 Philippe {Van Cappellen}, Y. Wang (1996) \n3.65e6 (mol l-1)-1 yr-1 Kazumi Ozaki, Shigeo Tajima, Eiichi Tajika (2011) \n54e6 (mol l-1)-1 yr-1 Stephen J. Romaniello, Louis A Derry (2010) \n\nParameters\n\nR_H2S_O2[Float64]=3650.0  ((mol m-3)-1 yr-1), default_value=3650.0, description=\"rate constant\"\nSIsotope[external, DataType]=PALEOboxes.ScalarData, default_value=PALEOboxes.ScalarData, allowed_values=Type[PALEOboxes.ScalarData, PALEOboxes.IsotopeLinear], description=\"disable / enable sulphur isotopes and specify isotope type\"\n\nMethods and Variables for default Parameters\n\ndo_redox_H2S_O2\nredox_H2S_O2  (mol O2 yr-1), VT_ReactProperty, description=\"oxygen consumption (+ve) by H2S oxidation\"\nO2_conc  (mol m-3), VT_ReactDependency, description=\"O2 concentration\"\nH2S_conc  (mol m-3), VT_ReactDependency, description=\"H2S concentration\"\nvolume  (m3), VT_ReactDependency, description=\"box fluid volume\"\nRateStoich_redox_H2S_O2\nredox_H2S_O2  (mol O2 yr-1), VT_ReactDependency, description=\"oxygen consumption (+ve) by H2S oxidation\"\n[O2_sms]  (), VT_ReactContributor, description=\"generated by RateStoich rate=redox_H2S_O2\"\n[H2S_sms]  (), VT_ReactContributor, description=\"generated by RateStoich rate=redox_H2S_O2\"\n[SO4_sms]  (), VT_ReactContributor, description=\"generated by RateStoich rate=redox_H2S_O2\"\n[TAlk_sms]  (), VT_ReactContributor, description=\"generated by RateStoich rate=redox_H2S_O2\"\ntotals\nredox_H2S_O2  (mol O2 yr-1), VT_ReactDependency, description=\"oxygen consumption (+ve) by H2S oxidation\"\nredox_H2S_O2_total  (mol O2 yr-1), VT_ReactProperty, description=\"total oxygen consumption (+ve) by H2S oxidation\"\n\n\n\n\n\n","category":"type"},{"location":"PALEOaqchem Reactions/#PALEOaqchem.SecondaryRedox.ReactionRedoxCH4_O2","page":"PALEOaqchem Reactions","title":"PALEOaqchem.SecondaryRedox.ReactionRedoxCH4_O2","text":"ReactionRedoxCH4_O2\n\nMethane oxidation by oxygen.\n\nRate R_CH4_O2 Units Ref Notes\n1e10 (mol l-1)-1 yr-1 Philippe {Van Cappellen}, Y. Wang (1996) \n10e6 (mol l-1)-1 yr-1 TODO \n\nParameters\n\nR_CH4_O2[Float64]=10000.0  ((mol m-3)-1 yr-1), default_value=10000.0, description=\"rate constant\"\nCIsotope[external, DataType]=PALEOboxes.ScalarData, default_value=PALEOboxes.ScalarData, allowed_values=Type[PALEOboxes.ScalarData, PALEOboxes.IsotopeLinear], description=\"disable / enable carbon isotopes and specify isotope type\"\n\nMethods and Variables for default Parameters\n\ndo_redox_CH4_O2\nredox_CH4_O2  (mol O2 yr-1), VT_ReactProperty, description=\"oxygen consumption (+ve) by CH4 oxidation\"\nO2_conc  (mol m-3), VT_ReactDependency, description=\"O2 concentration\"\nCH4_conc  (mol m-3), VT_ReactDependency, description=\"CH4 concentration\"\nvolume  (m3), VT_ReactDependency, description=\"box fluid volume\"\nRateStoich_redox_CH4_O2\nredox_CH4_O2  (mol O2 yr-1), VT_ReactDependency, description=\"oxygen consumption (+ve) by CH4 oxidation\"\n[O2_sms]  (), VT_ReactContributor, description=\"generated by RateStoich rate=redox_CH4_O2\"\n[CH4_sms]  (), VT_ReactContributor, description=\"generated by RateStoich rate=redox_CH4_O2\"\n[DIC_sms]  (), VT_ReactContributor, description=\"generated by RateStoich rate=redox_CH4_O2\"\ntotals\nredox_CH4_O2  (mol O2 yr-1), VT_ReactDependency, description=\"oxygen consumption (+ve) by CH4 oxidation\"\nredox_CH4_O2_total  (mol O2 yr-1), VT_ReactProperty, description=\"total oxygen consumption (+ve) by CH4 oxidation\"\n\n\n\n\n\n","category":"type"},{"location":"PALEOaqchem Reactions/#PALEOaqchem.SecondaryRedox.ReactionRedoxCH4_SO4","page":"PALEOaqchem Reactions","title":"PALEOaqchem.SecondaryRedox.ReactionRedoxCH4_SO4","text":"ReactionRedoxCH4_SO4\n\nMethane oxidation by sulphate (anaerobic methane oxidation).\n\nRate R_CH4_SO4 Units Ref Notes\n1e4 (mol l-1)-1 yr-1 Philippe {Van Cappellen}, Y. Wang (1996) k17\n\nParameters\n\nR_CH4_SO4[Float64]=10.0  ((mol m-3)-1 yr-1), default_value=10.0, description=\"rate constant\"\nCIsotope[external, DataType]=PALEOboxes.ScalarData, default_value=PALEOboxes.ScalarData, allowed_values=Type[PALEOboxes.ScalarData, PALEOboxes.IsotopeLinear], description=\"disable / enable carbon isotopes and specify isotope type\"\n\nMethods and Variables for default Parameters\n\ndo_redox_CH4_SO4\nredox_CH4_SO4  (mol CH4/SO4 yr-1), VT_ReactProperty, description=\"sulphate consumption (+ve) by CH4 oxidation\"\nSO4_conc  (mol m-3), VT_ReactDependency, description=\"SO4 concentration\"\nCH4_conc  (mol m-3), VT_ReactDependency, description=\"CH4 concentration\"\nvolume  (m3), VT_ReactDependency, description=\"box fluid volume\"\nRateStoich_redox_CH4_SO4\nredox_CH4_SO4  (mol CH4/SO4 yr-1), VT_ReactDependency, description=\"sulphate consumption (+ve) by CH4 oxidation\"\n[SO4_sms]  (), VT_ReactContributor, description=\"generated by RateStoich rate=redox_CH4_SO4\"\n[CH4_sms]  (), VT_ReactContributor, description=\"generated by RateStoich rate=redox_CH4_SO4\"\n[DIC_sms]  (), VT_ReactContributor, description=\"generated by RateStoich rate=redox_CH4_SO4\"\n[H2S_sms]  (), VT_ReactContributor, description=\"generated by RateStoich rate=redox_CH4_SO4\"\n[TAlk_sms]  (), VT_ReactContributor, description=\"generated by RateStoich rate=redox_CH4_SO4\"\ntotals\nredox_CH4_SO4  (mol CH4/SO4 yr-1), VT_ReactDependency, description=\"sulphate consumption (+ve) by CH4 oxidation\"\nredox_CH4_SO4_total  (mol CH4/SO4 yr-1), VT_ReactProperty, description=\"total sulphate consumption (+ve) by CH4 oxidation\"\n\n\n\n\n\n","category":"type"},{"location":"PALEOaqchem Reactions/#Carbonate-chemistry","page":"PALEOaqchem Reactions","title":"Carbonate chemistry","text":"","category":"section"},{"location":"PALEOaqchem Reactions/","page":"PALEOaqchem Reactions","title":"PALEOaqchem Reactions","text":"CarbChem.ReactionCO2SYS","category":"page"},{"location":"PALEOaqchem Reactions/#Isotope-systems","page":"PALEOaqchem Reactions","title":"Isotope systems","text":"","category":"section"},{"location":"PALEOaqchem Reactions/#Boron","page":"PALEOaqchem Reactions","title":"Boron","text":"","category":"section"},{"location":"PALEOaqchem Reactions/","page":"PALEOaqchem Reactions","title":"PALEOaqchem Reactions","text":"Boron.ReactionBoronIsotope","category":"page"},{"location":"PALEOaqchem Reactions/#PALEOaqchem.Boron.ReactionBoronIsotope","page":"PALEOaqchem Reactions","title":"PALEOaqchem.Boron.ReactionBoronIsotope","text":"ReactionBoronIsotope\n\nCalculate d11B for aqueous B(OH)4- and B(OH)3 species from mass balance, given  d11B for total B,  total B concentration B_conc and B(OH)4- concentration BOH4_conc.\n\nSee eg Richard E Zeebe, Dieter a. Wolf-Gladrow (2001), p220.\n\nParameters\n\nalphaB[Float64]=1.0272, default_value=1.0272, description=\"isotopic fractionation factor B(OH)4m <-> B(OH)3\"\nBIsotope[external, UnionAll]=PALEOboxes.IsotopeLinear, default_value=PALEOboxes.IsotopeLinear, allowed_values=Type[PALEOboxes.ScalarData, PALEOboxes.IsotopeLinear], description=\"disable / enable boron isotopes and specify isotope type\"\n\nMethods and Variables\n\ndo_boron_isotope\nB_conc  (mol m-3), VT_ReactDependency, description=\"B concentration\"\nBOH4_conc  (mol m-3), VT_ReactDependency, description=\"B(OH)4m concentration\"\nB_delta  (per mil), VT_ReactDependency, description=\"d11B delta for total B\"\nBOH4_delta  (per mil), VT_ReactProperty, description=\"d11B delta for B(OH)4- species\"\n\n\n\n\n\n","category":"type"},{"location":"#PALEOaqchem.jl","page":"PALEOaqchem.jl","title":"PALEOaqchem.jl","text":"","category":"section"},{"location":"","page":"PALEOaqchem.jl","title":"PALEOaqchem.jl","text":"Aquatic biogeochemistry for PALEO framework","category":"page"}]
}
