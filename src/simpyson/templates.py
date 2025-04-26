# Predefined pulse sequences

# No pulse sequence
no_pulse = """
proc pulseq {} {
    global par
    acq_block {
    delay $par(tsw)
    }
}
"""

# 90 degree pulse
pulse_90 = """
proc pulseq {} {
    global par
    acq_block {
    pulse $par(pH) $par(plH) $par(phH)  
    delay $par(tsw)
    }
}
"""