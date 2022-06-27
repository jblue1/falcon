import subprocess
import tensorflow as tf
import numpy as np
import matplotlib.pyplot as plt
from tensorflow import keras
import os
import sys
from numpy import pi
from scipy.stats import norm
from numpy.random import MT19937, RandomState, SeedSequence
import time

os.environ["CUDA_VISIBLE_DEVICES"] = ""
sys.path.append("../../models/cWGAN/Run_2021-07-30_0/")
import cWGAN

save_dir = "../../models/cWGAN/Run_2021-07-30_0/"
assert os.path.isdir(save_dir)


def print_bins(bin, i):
    print("bin {}".format(i))
    print(
        " Lower pPt: {:.3f} Upper pPt: {:.3f} Lower rPt: {:.3f} Upper rPt {:.3f}".format(
            bin[0], bin[8], bin[4], bin[12]
        )
    )
    print(
        " Lower pEta: {:.3f} Upper pEta: {:.3f} Lower rEta: {:.3f} Upper rEta {:.3f}".format(
            bin[1], bin[9], bin[5], bin[13]
        )
    )
    print(
        " Lower pPhi: {:.5f} Upper pPhi: {:.5f} Lower rPhi: {:.3f} Upper rPhi {:.3f}".format(
            bin[2], bin[10], bin[6], bin[14]
        )
    )
    print(
        " Lower pM: {:.3f} Upper pM: {:.3f} Lower rMi: {:.3f} Upper rM {:.3f}".format(
            bin[3], bin[11], bin[7], bin[15]
        )
    )
    print("")


def load_data(num_data_points=-1):
    np.random.seed(1)
    data = np.loadtxt(
        "../../data/processed/newPartonMatchedJetsNoRecoPtCutFixRapMassEventNumberJetTagged.txt",
        skiprows=2,
    )
    # print(data.shape)
    if num_data_points > 0:
        data = data[:num_data_points]
    np.random.shuffle(data)

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
    reco_data = (reco_data - reco_mean) / reco_std
    return (
        normalized_data,
        normalized_mean,
        normalized_std,
        data,
        reco_data,
        reco_mean,
        reco_std,
    )


def load_bins(normalized_mean, normalized_std, reco_mean, reco_std):
    #bins = np.loadtxt("../../data/processed/bins.txt", dtype=np.float32)
    bins = np.loadtxt("bins.txt")
    bins[:, :4] = (bins[:, :4] * normalized_std) + normalized_mean
    bins[:, 4:8] = (bins[:, 4:8] * reco_std) + reco_mean
    bins[:, 8:12] = (bins[:, 8:12] * normalized_std) + normalized_mean
    bins[:, 12:] = (bins[:, 12:] * reco_std) + reco_mean
    bins[:, 0] = 10 ** bins[:, 0]
    bins[:, 3] = 10 ** bins[:, 3]
    bins[:, 4] = 10 ** bins[:, 4]
    bins[:, 7] = 10 ** bins[:, 7]
    bins[:, 8] = 10 ** bins[:, 8]
    bins[:, 11] = 10 ** bins[:, 11]
    bins[:, 12] = 10 ** bins[:, 12]
    bins[:, 15] = 10 ** bins[:, 15]
    maxes = np.amax(bins[:, 8:], axis=0)
    mins = np.amin(bins[:, :8], axis=0)
    # make the mins and maxes for each dim +/- inf to ensure when we go back
    # and put everything back in the bins, we get all the data
    for i in range(bins.shape[0]):
        for j in range(bins.shape[1] // 2):
            if np.allclose(bins[i, j], mins[j]):
                bins[i, j] -= np.inf
            if np.allclose(bins[i, j + bins.shape[1] // 2], maxes[j]):
                bins[i, j + bins.shape[1] // 2] += np.inf
    return bins


def sort_data_into_bins(num_data_points, num_bins, data, bins, offset):
    true_bins_hist = np.zeros(num_bins)
    num_unbinned_jets = 0
    for i in range(num_data_points):
        found_bin = False
        for j in range(num_bins):
            if (
                data[offset * num_data_points + i, 0] >= bins[j, 0]
                and data[offset * num_data_points + i, 0] < bins[j, 8]
                and data[offset * num_data_points + i, 1] >= bins[j, 1]
                and data[offset * num_data_points + i, 1] < bins[j, 9]
                and data[offset * num_data_points + i, 2] >= bins[j, 2]
                and data[offset * num_data_points + i, 2] < bins[j, 10]
                and data[offset * num_data_points + i, 3] >= bins[j, 3]
                and data[offset * num_data_points + i, 3] < bins[j, 11]
                and data[offset * num_data_points + i, 4] >= bins[j, 4]
                and data[offset * num_data_points + i, 4] < bins[j, 12]
                and data[offset * num_data_points + i, 5] >= bins[j, 5]
                and data[offset * num_data_points + i, 5] < bins[j, 13]
                and data[offset * num_data_points + i, 6] >= bins[j, 6]
                and data[offset * num_data_points + i, 6] < bins[j, 14]
                and data[offset * num_data_points + i, 7] >= bins[j, 7]
                and data[offset * num_data_points + i, 7] < bins[j, 15]
            ):
                true_bins_hist[j] = true_bins_hist[j] + 1
                break
    print("Number of unbinned jets: ", num_unbinned_jets)
    return true_bins_hist


def main():
    rs = RandomState(MT19937(SeedSequence(123456789)))
    checkpoints = np.arange(0, 35000, 5000)
    task_id = int(sys.argv[1])
    print("Running task: ", task_id)
    num_tasks = int(sys.argv[2])
    data_file = "./data.txt"
    if task_id == 0:
        print("The task id was 0!")
        with open(data_file, 'w') as f:
            print("checkpoint num_bins NSIB_true_true NSIB_true_pred", file=f)
    else:
        while not os.path.exists(data_file):
            time.sleep(0.1)

    local_checkpoints = checkpoints[task_id::num_tasks]
    print("Task ", task_id, " will process checkpoints ", local_checkpoints)
    for checkpoint in local_checkpoints:
        print(f"Starting checkoint {checkpoint}")
        cwgan = cWGAN.cWGAN(10, "RMSprop", 0.000002, 0.00001, 10, False, "", 0)
        cwgan.generator.load_weights(save_dir + f"/training_checkpoints/gen_{str(checkpoint)}")
        (
            normalized_data,
            normalized_mean,
            normalized_std,
            data,
            reco_data,
            reco_mean,
            reco_std,
        ) = load_data(-1)
        predict = np.array(cwgan.make_generator_predictions(normalized_data))
        predict = predict * reco_std + reco_mean
        predict[:, 0] = 10 ** predict[:, 0]
        predict[:, 3] = 10 ** predict[:, 3]
        dataFile = "./shuffledData.txt"
        binsFile = "./bins.txt"
        thirds = len(data) // 3
        num_bins = 200
        print("Target points per bin: {}".format(thirds / num_bins))
        np.savetxt(dataFile, np.hstack((normalized_data, reco_data)))
        # if (kk==0):
        #    subprocess.run(["make", "makeKDTreeBins.out"])
        subprocess.run(
            ["../../bin/makeKDTreeBins.out", dataFile, binsFile, str(thirds), str(num_bins)]
        )

        bins = load_bins(normalized_mean, normalized_std, reco_mean, reco_std)

        num_data_points = thirds
        orig_bins_hist = sort_data_into_bins(num_data_points, num_bins, data, bins, 0)
        print("orig_bins_hist = ", orig_bins_hist)
        true_bins_hist = sort_data_into_bins(num_data_points, num_bins, data, bins, 1)
        print("true_bins_hist = ", true_bins_hist)
        Np = num_data_points
        Nq = num_data_points
        Nr = num_data_points
        sumIBp = num_data_points / num_bins
        Pp = sumIBp / num_data_points

        p_vals = []

        for i in range(num_bins):
            sumIBq = true_bins_hist[i]
            Pq = sumIBq / num_data_points
            P = (sumIBp + sumIBq) / (2 * num_data_points)
            Se = np.sqrt(P * (1 - P) * ((1 / Np) + (1 / Nq)))
            z_score = (Pp - Pq) / Se
            p_vals.append(2 * (1 - norm.cdf(np.abs(z_score))))

        s = 0
        for i in range(len(p_vals)):
            if p_vals[i] < 0.05:
                s += 1
        
        nsib_true_true = s * 100 / len(p_vals)
        print(
            "Proportion of statistically distinct bins: {:.3f}".format(
                s * 100 / len(p_vals)
            )
        )
        #propDisBins1[kk] = s * 100 / len(p_vals)

        p_vals.clear()

        pred_bins_hist = np.zeros(num_bins)
        for i in range(num_data_points):
            for j in range(num_bins):
                if (
                    data[2 * num_data_points + i, 0] >= bins[j, 0]
                    and data[2 * num_data_points + i, 0] < bins[j, 8]
                    and data[2 * num_data_points + i, 1] >= bins[j, 1]
                    and data[2 * num_data_points + i, 1] < bins[j, 9]
                    and data[2 * num_data_points + i, 2] >= bins[j, 2]
                    and data[2 * num_data_points + i, 2] < bins[j, 10]
                    and data[2 * num_data_points + i, 3] >= bins[j, 3]
                    and data[2 * num_data_points + i, 3] < bins[j, 11]
                    and predict[2 * num_data_points + i, 0] >= bins[j, 4]
                    and predict[2 * num_data_points + i, 0] < bins[j, 12]
                    and predict[2 * num_data_points + i, 1] >= bins[j, 5]
                    and predict[2 * num_data_points + i, 1] < bins[j, 13]
                    and predict[2 * num_data_points + i, 2] >= bins[j, 6]
                    and predict[2 * num_data_points + i, 2] < bins[j, 14]
                    and predict[2 * num_data_points + i, 3] >= bins[j, 7]
                    and predict[2 * num_data_points + i, 3] < bins[j, 15]
                ):
                    pred_bins_hist[j] = pred_bins_hist[j] + 1

        print("pred_bins_hist = ", pred_bins_hist)
        for i in range(num_bins):
            sumIBq = true_bins_hist[i]
            sumIBr = pred_bins_hist[i]
            Pq = sumIBq / num_data_points
            Pr = sumIBr / num_data_points
            P = (sumIBr + sumIBq) / (2 * num_data_points)
            Se = np.sqrt(P * (1 - P) * ((1 / Nr) + (1 / Nq)))
            z_score = (Pr - Pq) / Se
            p_vals.append(2 * (1 - norm.cdf(np.abs(z_score))))

        s = 0
        for i in range(len(p_vals)):   
            if p_vals[i] < 0.05:
                s += 1
        nsib_true_pred = s * 100 / len(p_vals)
        print(
            "Proportion of statistically distinct bins: {:.3f}".format(
                s * 100 / len(p_vals)
            )
        )
        print("Should be adding line to file")
        with open(data_file, 'a') as f:
            print(f"{checkpoint} {num_bins} {nsib_true_true} {nsib_true_pred}", file=f)
    #propDisBins2[kk] = s * 100 / len(p_vals)

    """
    if kk == 5:
        plt.rc("font", family="serif")
        plt.rc("text", usetex=True)
        fig = plt.figure(figsize=(8, 6), dpi=100)
        ax = fig.add_subplot(111)
        ax.set_title("Proportion of Samples in Bins formed from KDTrees")
        ax.bar(
            np.arange(len(true_bins_hist)),
            true_bins_hist / num_data_points,
            width=1,
            label="Real Data",
        )
        ax.bar(
            np.arange(len(true_bins_hist)),
            pred_bins_hist / num_data_points,
            width=1,
            alpha=0.5,
            label="Pred. Data",
        )
        ax.set_xlabel("Bin Number")
        ax.set_ylabel("Proportion of Samples")
        ax.legend()
        fig.savefig("./fullHist1.png")
    """
    # np.savetxt("./propDisBins1.txt", propDisBins1)
    # np.savetxt("./propDisBins2.txt", propDisBins2)


if __name__ == "__main__":
    main()
