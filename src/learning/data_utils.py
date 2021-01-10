def load_jet_data(data_path, batch_size):
    """
    Load and normalize the parton data.
    data_path - path to txt file with jet 4-momenta
    returns - tuple of numpy arrays (parton_data, reco_data)
    """
    data = np.loadtxt(data_path, skiprows=2)
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


def one_hot_encode(number):
    """
    Return one-hot-encoded vector to serve as a label for MNIST data, for
    example 
    one_hot_encode(4)
    would return
    [0, 0, 0, 0, 1, 0, 0, 0, 0, 0]
    number - digit (must be 0 through 9)
    returns - one hot encoded vector
    """
    if (number < 0 or number > 9):
        sys.exit(1)
    label = np.zeros(10)
    label[number] = 1
    return label


def load_mnist_data():
    """
    Load and normalize MNIST data
    returns - tuple of numpy arrays (labels, images)
    """
    (train_images, train_labels), (_, _) = tf.keras.datasets.mnist.load_data()
    train_images = train_images.reshape(train_images.shape[0], 28, 28, 1).astype('float32')
    train_images = (train_images - 127.5) / 127.5 # Normalize the images to [-1, 1]

    one_hot_labels = np.zeros((train_labels.shape[0], 10))
    for i in range(len(train_labels)):
        one_hot_labels[i, :] = one_hot_encode(train_labels[i])

    return (train_labels, one_hot_labels, train_images)
