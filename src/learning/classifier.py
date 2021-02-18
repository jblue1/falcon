import tensorflow as tf
from tensorflow import keras
import data_utils
import file_utils
import os


def build_model():
    model = keras.Sequential()
    model.add(keras.Input(shape=(2,)))
    model.add(keras.layers.Dense(16, activation="relu"))
    model.add(keras.layers.Dense(16, activation="relu"))
    model.add(keras.layers.Dense(1, activation="sigmoid"))
    return model

class Trainer:
    def __init__(self, params_dict):
        """Constructor

        Args:
        """
        self.save_dir = file_utils.make_save_directory("Classifier")
        self.model = build_model()
        self.model.compile(
            optimizer=tf.keras.optimizers.Adam(learning_rate=params_dict['lr']),
            loss=tf.keras.losses.BinaryCrossentropy(),
            metrics=[tf.keras.metrics.BinaryAccuracy()],
        )

        self.batch_size = params_dict["batch_size"]
        self.data = data_utils.load_classifier_data(params_dict["data_path"])
        self.checkpoint_path = os.path.join(self.save_dir, "training/cp.cpkt")
        self.epochs = params_dict["epochs"]

    def train(self):
        cp_callback = tf.keras.callbacks.ModelCheckpoint(
            filepath=self.checkpoint_path,
            save_weights_only=True,
            verbose=1,
            save_best_only=True,
        )

        history = self.model.fit(
            self.data[:, 0:2],
            self.data[:, 2],
            batch_size=self.batch_size,
            epochs=self.epochs,
            validation_split=0.2,
            callbacks=[cp_callback],
        )

        return history

    def save_losses(self, history):
        training_loss = history.history["loss"]
        val_loss = history.history["val_loss"]
        loss_dict = {"Training Loss": training_loss, "Validation Loss": val_loss}
        file_utils.save_losses(self.save_dir, loss_dict)

    def save_model(self):
        file_utils.save_network(self.save_dir, model_path="./classifier.py")

    def save_params(self, params_dict):
        file_utils.save_params(self.save_dir, params_dict)

    