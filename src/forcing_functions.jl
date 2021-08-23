function convert_forcing_thetal(param_set, p0, qt, qv, T, qt_tendency, T_tendency)
    # TODO - should pp also be passed here?
    return T_tendency / TD.exner_given_pressure(param_set, p0)
end
