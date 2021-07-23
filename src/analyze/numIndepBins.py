import tensorflow as tf
import numpy as np
import matplotlib.pyplot as plt
from tensorflow import keras
import os
import sys
from numpy import pi

os.environ["CUDA_VISIBLE_DEVICES"] = ""
sys.path.append("./models/cWGAN/Run_2021-07-15_1/")
import cWGAN

save_dir = "./models/cWGAN/Run_2021-07-15_1/"
assert os.path.isdir(save_dir)


def load_data():
    data = np.loadtxt(
        "./data/processed/newPartonMatchedJetsNoRecoPtCutFixRapMass.txt", skiprows=2
    )

    normalized_data = np.zeros((len(data), 4))
    reco_data = np.zeros((len(data), 4))

    # for log scaling
    np.log10(data[:, 0], out=normalized_data[:, 0])
    np.log10(data[:, 3], out=normalized_data[:, 3])
    np.log10(data[:, 4], out=reco_data[:, 0])
    np.log10(data[:, 7], out=reco_data[:, 3])

    normalized_data[:, 1:3] = data[:, 1:3]
    reco_data[:, 1:3] = data[:, 5:7]

    normalized_mean = np.mean(normalized_data, axis=0)
    normalized_std = np.std(normalized_data, axis=0)
    reco_mean = np.mean(reco_data, axis=0)
    reco_std = np.std(reco_data, axis=0)

    normalized_data = (normalized_data - normalized_mean) / normalized_std
    return normalized_data, data, reco_mean, reco_std


def main():
    normalized_data, data, reco_mean, reco_std = load_data()
    cwgan = cWGAN.cWGAN(10, "RMSprop", 0.000002, 0.00001, 10, False, "", 0)
    cwgan.generator.load_weights(save_dir + '/training_checkpoints/gen_175000')
    predict = np.array(cwgan.make_generator_predictions(normalized_data))
    predict = predict * reco_std + reco_mean
    predict[:, 0] = 10**predict[:, 0]
    predict[:, 3] = 10**predict[:, 3]
    bins = np.loadtxt('./data/processed/testbins.txt', dtype=np.float32)
    maxes = np.amax(bins[:, 8:], axis=0)
    mins = np.amin(bins[:, :8], axis=0)
    for i in range(bins.shape[0]):
        for j in range(bins.shape[1]//2):
            if np.allclose(bins[i, j], mins[j]):
                bins[i, j] -= 1
            if np.allclose(bins[i, j+bins.shape[1]//2], maxes[j]):
                bins[i, j+bins.shape[1]//2] += 1
    
    num_bins = len(bins)
    num_data_points = 10000
    bins_hist = np.zeros(num_bins)
    Np = num_data_points / num_bins
    Pp = Np / num_data_points
    print("Np: {}".format(Np))
    print("Pp: {}".format(Pp))
    
    for i in range(num_data_points):
        #if (i % (num_data_points/10) == 0):
         #   print(i/10)
        for j in range(num_bins):
            if (data[num_data_points+i, 0] >= bins [j, 0] and data[num_data_points+i, 0] < bins[j, 8] and
                data[num_data_points+i, 1] >= bins [j, 1] and data[num_data_points+i, 1] < bins[j, 9] and
                data[num_data_points+i, 2] >= bins [j, 2] and data[num_data_points+i, 2] < bins[j, 10] and
                data[num_data_points+i, 3] >= bins [j, 3] and data[num_data_points+i, 3] < bins[j, 11] and
                data[num_data_points+i, 4] >= bins [j, 4] and data[num_data_points+i, 4] < bins [j, 12] and
                data[num_data_points+i, 5] >= bins [j, 5] and data[num_data_points+i, 5] < bins [j, 13] and
                data[num_data_points+i, 6] >= bins [j, 6] and data[num_data_points+i, 6] < bins [j, 14] and
                data[num_data_points+i, 7] >= bins [j, 7] and data[num_data_points+i, 7] < bins [j, 15]):
                bins_hist[j] = bins_hist[j] + 1


    for i in range(num_bins):
        Nq = bins_hist[i]
        Pq = Nq / num_data_points
        print("Nq: {}".format(Nq))
        print("Pq: {}".format(Pq))
        P = (Np + Nq) / (2*num_data_points)
        print("P: {}".format(P))
        Se = np.sqrt(P*(1-P)*( (1/Np)+(1/Nq) ))
        print("Se: {}".format(Se))
        z_score = (Pp - Pq / Se)
        print(z_score)


    print(bins_hist)
if __name__ == "__main__":
    main()
