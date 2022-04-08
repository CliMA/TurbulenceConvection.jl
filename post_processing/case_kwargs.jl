import OrderedCollections
case_kwargs = OrderedCollections.OrderedDict()

#! format: off
case_kwargs["Bomex"] = (;
    plot_comparison = true,
    t_start = 4 * 3600,
    t_stop = 6 * 3600,
)

case_kwargs["life_cycle_Tan2018"] = (;
    plot_comparison = true,
    t_start = 4 * 3600,
    t_stop = 6 * 3600)

case_kwargs["Soares"] = (;
    plot_comparison = true,
    t_start = 6 * 3600,
    t_stop = 8 * 3600,)

case_kwargs["Rico"] = (;
    plot_comparison = true,
    t_start = 22 * 3600,
    t_stop = 24 * 3600)

case_kwargs["Nieuwstadt"] = (;
    plot_comparison = true,
    t_start = 6 * 3600,
    t_stop = 8 * 3600)

case_kwargs["TRMM_LBA"] = (;
    plot_comparison = true,
    t_start = 4 * 3600,
    t_stop = 6 * 3600)

case_kwargs["ARM_SGP"] = (;
    plot_comparison = true,
    t_start = 8 * 3600,
    t_stop = 11 * 3600)

case_kwargs["GATE_III"] = (;
    plot_comparison = true,
    t_start = 4 * 3600,
    t_stop = 6 * 3600,
    skip_comparison = true # TODO: remove this option
)

case_kwargs["DYCOMS_RF01"] = (;
    plot_comparison = true,
    t_start = 2 * 3600,
    t_stop = 4 * 3600)

case_kwargs["DYCOMS_RF02"] = (;
    plot_comparison = true,
    t_start = 4 * 3600,
    t_stop = 6 * 3600)

case_kwargs["GABLS"] = (;
    plot_comparison = true,
    t_start = 7 * 3600,
    t_stop = 9 * 3600)

case_kwargs["SP"] = (;
    plot_comparison = true,
    t_start = 0,
    t_stop = 2 * 3600)

case_kwargs["DryBubble"] = (;
    plot_comparison = true,
    t_start = 900,
    t_stop = 1000)

case_kwargs["LES_driven_SCM"] = (;
    ds_les_filename = joinpath(NameList.les_driven_scm_data_folder(), "Stats.cfsite23_HadGEM2-A_amip_2004-2008.07.nc"),
    plot_comparison = true,
    t_start = 3 * 3600,
    t_stop = 6 * 3600,
)

#! format: on
