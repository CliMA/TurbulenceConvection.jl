function set_bcs(self::EnvironmentVariable, Gr::Grid)
    start_low = Gr.gw
    start_high = Gr.nzg - Gr.gw

    if self.name == "w"
        self.values[start_high] = 0.0
        self.values[start_low] = 0.0
        @inbounds for kk in xrange(1, Gr.gw)
            k = kk - 1
            self.values[start_high + k] = -self.values[start_high - k]
            self.values[start_low - k] = -self.values[start_low + k]
        end
    else
        @inbounds for kk in xrange(Gr.gw)
            k = kk - 1
            self.values[start_high + k + 1] = self.values[start_high - k]
            self.values[start_low - k] = self.values[start_low + 1 + k]
        end
    end
end

function set_bcs(self::EnvironmentVariable_2m, Gr::Grid)
    start_low = Gr.gw
    start_high = Gr.nzg - Gr.gw

    @inbounds for kk in xrange(Gr.gw)
        k = kk - 1
        self.values[start_high + k + 1] = self.values[start_high - k]
        self.values[start_low - k] = self.values[start_low + 1 + k]
    end
end

function set_bcs(self::RainVariable, Gr::Grid)
    @inbounds for k in xrange(Gr.gw)
        self.values[Gr.nzg - Gr.gw + k] = self.values[Gr.nzg - Gr.gw - 1 - k]
        self.values[Gr.gw - 1 - k] = self.values[Gr.gw + k]
    end
    return
end

function set_bcs(self::VariablePrognostic, Gr::Grid)
    start_low = Gr.gw
    start_high = Gr.nzg - Gr.gw

    if self.bc == "sym"
        @inbounds for kk in xrange(Gr.gw)
            k = kk - 1
            self.values[start_high + k + 1] = self.values[start_high - k]
            self.values[start_low - k] = self.values[start_low + 1 + k]

            self.mf_update[start_high + k + 1] = self.mf_update[start_high - k]
            self.mf_update[start_low - k] = self.mf_update[start_low + 1 + k]

            self.new[start_high + k + 1] = self.new[start_high - k]
            self.new[start_low - k] = self.new[start_low + 1 + k]
        end
    else
        self.values[start_high] = 0.0
        self.values[start_low] = 0.0

        self.mf_update[start_high] = 0.0
        self.mf_update[start_low] = 0.0

        self.new[start_high] = 0.0
        self.new[start_low] = 0.0

        @inbounds for kk in xrange(1, Gr.gw)
            k = kk - 1
            self.values[start_high + k] = -self.values[start_high - k]
            self.values[start_low - k] = -self.values[start_low + k]

            self.mf_update[start_high + k] = -self.mf_update[start_high - k]
            self.mf_update[start_low - k] = -self.mf_update[start_low + k]

            self.new[start_high + k] = -self.new[start_high - k]
            self.new[start_low - k] = -self.new[start_low + k]
        end
    end

    return
end

function set_bcs(self::VariableDiagnostic, Gr::Grid)
    start_low = Gr.gw
    start_high = Gr.nzg - Gr.gw + 1

    if self.bc == "sym"
        @inbounds for kk in xrange(Gr.gw)
            k = kk - 1
            self.values[start_high + k] = self.values[start_high - 1]
            self.values[start_low - k] = self.values[start_low + 1]
        end

    else
        self.values[start_high] = 0.0
        self.values[start_low] = 0.0
        @inbounds for kk in xrange(1, Gr.gw)
            k = kk - 1
            self.values[start_high + k] = 0.0  #-self.values[start_high - k ]
            self.values[start_low - k] = 0.0 #-self.values[start_low + k ]
        end
    end

    return
end

function set_bcs(self::UpdraftVariable, Gr::Grid)
    start_low = Gr.gw
    start_high = Gr.nzg - Gr.gw

    n_updrafts = size(self.values)[1]

    if self.name == "w"
        @inbounds for i in xrange(n_updrafts)
            self.values[i, start_high] = 0.0
            self.values[i, start_low] = 0.0
            @inbounds for kk in xrange(1, Gr.gw)
                k = kk - 1
                self.values[i, start_high + k] = -self.values[i, start_high - k]
                self.values[i, start_low - k] = -self.values[i, start_low + k]
            end
        end
    else
        @inbounds for kk in xrange(Gr.gw)
            k = kk - 1
            @inbounds for i in xrange(n_updrafts)
                self.values[i, start_high + k + 1] = self.values[i, start_high - k]
                self.values[i, start_low - k] = self.values[i, start_low + 1 + k]
            end
        end
    end
    return
end
