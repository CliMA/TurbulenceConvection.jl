function convert_forcing_thetal(p0, qt, qv, T, qt_tendency, T_tendency)
    return T_tendency / exner_c(p0)
end
