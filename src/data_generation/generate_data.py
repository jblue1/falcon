"""
Script that generates a data set of parton jet and reco jet pairs. This involves
generating events with pythia (handled by event_gen.cc) and simulating a detector with
delphes (handled by process_delphes_output.cc). This script assumes that you have compiled
both of those files and that the executables live in the bin/ directory. You also need
to manually create a subdirectory
data/<date>
where the date is in the form
<year>-<month>-<day>
Additionally, you need to have a file called pythia_settings.txt inside this directory
with the settings for the pythia event (an example is included in the repo).

The output of this script are the following files:
partonJets.txt - a list of parton jet 4 momenta along with the event number
genJets.txt - a list of gen jet 4 momenta along with the event number
delphesOut.root - the output of the Delphes ExRootTreeWriter
matchedJets.txt - a list of matched parton and reco jet 4 momenta
pythia_settings.txt - the pythia settings file will get copied over and the random seed
used by the pythia object will get appended as another line

These files will all be placed in the directory
data/<date>/<task-number>
"""

import os
import datetime
import sys
import subprocess
from random import randint


def main():
    if os.environ["ROOT_INCLUDE_PATH"] == "":
        sys.exit("Must set ROOT_INCLUDE_PATH to path/to/delphes/external")
    num_events = sys.argv[1]
    delphes_path = sys.argv[2]
    today = str(datetime.date.today())
    save_dir = "../../data/" + str(today) + "/"
    print("Checking if directory ", save_dir, " exists")
    if not os.path.exists(save_dir):
        sys.exit("Data directory does not exist")
    seed = randint(1, 900000000)

    # EVENT GENERATION
    command1 = (
        "../../bin/event_generation.out"
        + " "
        + str(num_events)
        + " "
        + save_dir
        + " "
        + delphes_path
        + " "
        + str(seed),
    )
    output1 = subprocess.run(
        command1,
        capture_output=True,
        shell=True,
    )
    print(output1.stdout.decode("utf-8"))
    print(output1.stderr.decode("utf-8"))
    if output1.returncode != 0:
        sys.exit("Event generation failed.")

    # PROCESS JETS
    command2 = "../../bin/process_delphes_output.out " + save_dir
    output2 = subprocess.run(
        command2,
        capture_output=True,
        shell=True,
    )
    print(output2.stdout.decode("utf-8"))
    print(output2.stderr.decode("utf-8"))
    if output2.returncode != 0:
        sys.exit("Jet processing failed.")


if __name__ == "__main__":
    main()