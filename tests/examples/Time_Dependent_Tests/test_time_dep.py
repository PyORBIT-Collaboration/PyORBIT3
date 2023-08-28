import os
import logging

from orbit.time_dep import time_dep
from orbit.time_dep.waveform import LinearMagnetWaveform
from orbit.core.bunch import Bunch


# this function is just to replace the DEBUG message that logger prints so that we can compare it to the python2.7 output
def trim_logoutput(astring):
    lines = astring.split("\n")
    modified_lines = [line.replace("DEBUG    orbit.time_dep.time_dep:time_dep.py:121 ", "") for line in lines]
    modified_string = "\n".join(modified_lines)
    return modified_string


def test_time_dep(caplog):
    caplog.set_level(logging.DEBUG)
    script_dir = os.path.dirname(__file__)
    rcs_1124 = os.path.join(script_dir, "rcs1124.dat")

    b = Bunch()
    b.addParticle(1.0e-3, 0.0, 0.0, 0.0, 0.0, 0.0)
    b.addParticle(0.0, 1.0e-3, 0.0, 0.0, 0.0, 0.0)
    b.addParticle(0.0, 0.0, 1.0e-3, 0.0, 0.0, 0.0)
    b.addParticle(0.0, 0.0, 0.0, 1.0e-3, 0.0, 0.0)
    b.addParticle(0.0, 0.0, 0.0, 0.0, 1.0e-3, 0.0)
    b.addParticle(0.0, 0.0, 0.0, 0.0, 0.0, 1.0e-3)
    b.compress()
    syncPart = b.getSyncParticle()
    energy = 1.0
    syncPart.kinEnergy(energy)

    time_dep_latt = time_dep.TIME_DEP_Lattice()
    time_dep_latt.readMAD(rcs_1124, "RCS")
    time_dep_latt.setLatticeOrder()

    WaveForm01 = LinearMagnetWaveform()
    WaveForm01.initialize(0, 1.0, 1.0, 0.97)

    time_dep_latt.setTimeDepNode("R1SF02_1", WaveForm01)
    time_dep_latt.setTurns(10)
    time_dep_latt.trackBunchTurns(b)

    expected_output = """trackBunchTurns time 8.68842277975e-07 in 0 turn
trackBunchTurns time 1.73768455595e-06 in 1 turn
trackBunchTurns time 2.60652683393e-06 in 2 turn
trackBunchTurns time 3.4753691119e-06 in 3 turn
trackBunchTurns time 4.34421138988e-06 in 4 turn
trackBunchTurns time 5.21305366785e-06 in 5 turn
trackBunchTurns time 6.08189594583e-06 in 6 turn
trackBunchTurns time 6.9507382238e-06 in 7 turn
trackBunchTurns time 7.81958050178e-06 in 8 turn\n"""

    log_output = trim_logoutput(caplog.text)
    assert log_output == expected_output
