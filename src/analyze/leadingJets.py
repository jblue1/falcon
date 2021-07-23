"""A script to analyze model predictions in an event-by-event fashion.
"""

import numpy as np
import matplotlib.pyplot as plt
import tensorflow as tf
from tensorflow import keras
import os
import sys

# ensure no GPU usage
os.environ["CUDA_VISIBLE_DEVICES"] = ""

# choose the model to use
sys.path.append("./models/cWGAN/Run_2021-07-15_1/")
import cWGAN

save_dir = "./models/cWGAN/Run_2021-07-15_1/"
assert os.path.isdir(save_dir)

def loadData():
    data = np.loadtxt(
        "./data/processed/newPartonMatchedJetsNoRecoPtCutFixRapMassEventNumber.txt",
        skiprows=2
    )

    np.log10(data[:, 0], out=data[:, 0])
    np.log10(data[:, 3], out=data[:, 3])
    np.log10(data[:, 4], out=data[:, 4])
    np.log10(data[:, 7], out=data[:, 7])

    normalized_mean = np.mean(data[:, :4], axis=0)
    normalized_std = np.std(data[:, :4], axis=0)
    reco_mean = np.mean(data[:, 4:8], axis=0)
    reco_std = np.std(data[:, 4:8], axis=0)

    data[:, :4] = (data[:, :4] - normalized_mean) / normalized_std
    data[:, 4] = 10**data[:, 4]
    data[:, 7] = 10**data[:, 7]
    
    return data, reco_mean, reco_std


def loadEvents(data):
    eventIndex = 0
    jetIndex = 0
    events = []
    jets = []
    for jetIndex in range(len(data)):
        if abs(data[jetIndex, 8] - eventIndex) < 1e-6:
            jets.append(data[jetIndex, :].tolist())
        else:
            arr = np.asarray(jets.copy())
            arr.sort(axis=0)
            arr = np.flip(arr, 0)
            events.append(arr)
            jets.clear()
            jets.append(data[jetIndex, :].tolist())
            eventIndex += 1

    arr = np.asarray(jets.copy())
    arr.sort(axis=0)
    arr = np.flip(arr, 0)
    events.append(arr)

    return events


def main():
    cwgan = cWGAN.cWGAN(10, "RMSprop", 0.000002, 0.00001, 10, False, "", 0)
    cwgan.generator.load_weights(save_dir + '/training_checkpoints/gen_195000')
    data, reco_mean, reco_std = loadData()
    events = loadEvents(data)
    leadingJets = []
    secondLeadingJets = []
    predictedLeadingJets = []
    predictedSecondLeadingJets = []
    for event in events:
        if len(event) >= 6:
            inputs = np.vstack((event[0][:4], event[5][:4]))
            predict = np.array(cwgan.make_generator_predictions(inputs))
            predict = predict * reco_std + reco_mean
            predict[:, 0] = 10**predict[:, 0]
            predict[:, 3] = 10**predict[:, 3]
            predictedLeadingJets.append(predict[0, 0])
            predictedLeadingJets.append(predict[1, 0])
            #predictedSecondLeadingJets.append(predict[1, 0])
            leadingJets.append(event[0][4])
            leadingJets.append(event[5][4])
            #secondLeadingJets.append(event[5][4])
    print(len(leadingJets))
    bins = np.linspace(0, 400, 100)
    fig = plt.figure(figsize=(8,6), dpi=100)
    ax = fig.add_subplot(111)
    ax.hist(leadingJets, bins=bins, label="True")
    ax.hist(predictedLeadingJets, bins=bins, alpha=0.5, label="Predicted")
    ax.set_title(r"$p_T$ Distribution of First and Sixth leading $p_T$ Jet in each event")
    ax.legend()
    plt.savefig("./hist2.png")


if __name__ == "__main__":
    main()
