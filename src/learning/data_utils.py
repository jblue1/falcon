import numpy as np
import tensorflow as tf
import sys


def load_jet_data(data_path):
    """Load and normalize the jet data

    For each 4-momentum, scale Pt and E by subtracting the min and
    dividing by the max (i.e scale to [0, 1]), and scale eta and phi
    by subtracting by the mean and dividing by the std deviation 
    (i.e scale to mean of 0 and unit variance)
    Args:
        data_path (path-like): path to txt file with jet 4-momenta

    Returns:
        tuple: tuple of ndarrays (parton_data, reco_data)
    """
    data = np.loadtxt(data_path, skiprows=2, dtype=np.float32)
    partonPtMax = np.max(data[:, 0], axis=0)
    partonPtMin = np.min(data[:, 0], axis=0)
    partonMean = np.mean(data[:, 1:3], axis=0)
    partonStd = np.std(data[:, 1:3], axis=0)
    partonEMax = np.max(data[:, 3], axis=0)
    partonEMin = np.min(data[:, 3], axis=0)

    pfPtMax = np.max(data[:, 4], axis=0)
    pfPtMin = np.min(data[:, 4], axis=0)
    pfMean = np.mean(data[:, 5:7], axis=0)
    pfStd = np.std(data[:, 5:7], axis=0)
    pfEMax = np.max(data[:, 7], axis=0)
    pfEMin = np.min(data[:, 7], axis=0)

    data[:, 0] = (data[:, 0] - partonPtMin) / partonPtMax
    data[:, 1:3] = (data[:, 1:3] - partonMean) / partonStd
    data[:, 3] = (data[:, 3] - partonEMin) / partonEMax
    data[:, 4] = (data[:, 4] - pfPtMin) / pfPtMax
    data[:, 5:7] = (data[:, 5:7] - pfMean) / pfStd
    data[:, 7] = (data[:, 7] - pfEMin) / pfEMax

    np.random.shuffle(data)
    parton_data = data[:, :4]
    reco_data = data[:, 4:]
    return (parton_data, reco_data)


def load_jet_data_inverse_scaling(data_path):
    """Load and normalize the jet data

    For Pt and E, instead of scaling to [0, 1] by subtracting the min and adding the max,
    this function scales each by 1/Pt (or 1/E). 
    Args:
        data_path (path-like): path to txt file with jet 4-momenta

    Returns:
        tuple: tuple of ndarrays (parton_data, reco_data)
    """
    data = np.loadtxt(data_path, skiprows=2, dtype=np.float32)
    partonMean = np.mean(data[:, 1:3], axis=0)
    partonStd = np.std(data[:, 1:3], axis=0)

    pfMean = np.mean(data[:, 5:7], axis=0)
    pfStd = np.std(data[:, 5:7], axis=0)

    data[:, 0] = 1/data[:, 0]
    data[:, 1:3] = (data[:, 1:3] - partonMean) / partonStd
    data[:, 3] = 1/data[:, 3]
    data[:, 4] = 1/data[:, 4]
    data[:, 5:7] = (data[:, 5:7] - pfMean) / pfStd
    data[:, 7] = 1/data[:, 7]

    np.random.shuffle(data)
    parton_data = data[:, :4]
    reco_data = data[:, 4:]
    return (parton_data, reco_data)


def one_hot_encode(number):
    """Return a one-hot-encoded vector to serve as a label for MNIST data.

    For example, one_hot_encode(4) would return

    [0, 0, 0, 0, 1, 0, 0, 0, 0, 0]

    Args:
        number (int): the digit to encode

    Returns:
        ndarray: the one-hot-encoded vector
    """
    if number < 0 or number > 9:
        sys.exit(1)
    label = np.zeros(10)
    label[number] = 1.0
    return label


def load_mnist_data():
    """Load and normalize MNIST data

    Returns:
        tuple: tuple of ndarrays (labels, images)
    """
    (train_images, train_labels), (_, _) = tf.keras.datasets.mnist.load_data()
    train_images = train_images.reshape(train_images.shape[0], 28, 28, 1).astype(
        "float32"
    )
    train_images = (train_images - 127.5) / 127.5  # Normalize the images to [-1, 1]

    one_hot_labels = np.zeros((train_labels.shape[0], 10), dtype=np.float32)
    for i in range(len(train_labels)):
        one_hot_labels[i, :] = one_hot_encode(train_labels[i])

    return (one_hot_labels, train_images)


def concatenate_images_labels(images, labels):
    """Concatenate labels to images depthwise

    For example, say you are training on MNIST, so the images are
    BATCH_SIZEx28x28x1, and the labels are BATCH_SIZEx10x1. The output would be
    BATCH_SIZE x 28 x28 x 11

    Args:
        images (tf.Tensor): Batch of images
        labels (tf.Tensor): Batch of labels

    Returns:
        tf.Tensor: Batch with labels concatenated with the images
    """
    labels = tf.expand_dims(labels, axis=1)
    labels = tf.expand_dims(labels, axis=1)
    ones = tf.ones(images.shape)
    labels_for_concat = labels * ones
    return tf.concat((images, labels_for_concat), axis=3)


def test_concatenate_images_labels():
    """Test the concatenate_images_labels function with an example thats easy to solve by hand."""
    print("Testing concatenate labels")
    images = tf.constant([[[[2, 3], [4, 5]], [[6, 7], [8, 9]]]], dtype=tf.int8)

    images = tf.reshape(images, (2, 2, 2, 1))

    labels = tf.constant([[1, 0], [0, 1]], dtype=tf.int8)
    desired_output = tf.constant(
        [
            [[[2, 1, 0], [3, 1, 0]], [[4, 1, 0], [5, 1, 0]]],
            [[[6, 0, 1], [7, 0, 1]], [[8, 0, 1], [9, 0, 1]]],
        ],
        dtype=tf.int8,
    )

    images_and_labels = concatenate_images_labels(images, labels)
    tf.debugging.assert_equal(images_and_labels, desired_output)


def main():
    test_concatenate_images_labels()


if __name__ == "__main__":
    main()
