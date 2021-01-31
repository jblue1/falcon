import tensorflow as tf
from tensorflow import keras
import file_utils
import data_utils
import os


def make_model():
    model = keras.Sequential()
    model.add(keras.Input(shape=(4,)))
    model.add(keras.layers.Dense(128, activation="relu"))
    model.add(keras.layers.Dense(128, activation="relu"))
    model.add(keras.layers.Dense(128, activation="relu"))
    model.add(keras.layers.Dense(128, activation="relu"))
    model.add(keras.layers.Dense(128, activation="relu"))
    model.add(keras.layers.Dense(4))
    return model


class Trainer:
    def __init__(self, params_dict):
    #def __init__(self, lr, epochs, data_path, batch_size):
        self.save_dir = file_utils.make_save_directory("FCNN")
        self.model = make_model()
        self.model.compile(
            optimizer=tf.keras.optimizers.Adam(learning_rate=params_dict["lr"]),
            loss=tf.keras.losses.MeanAbsoluteError(),
            metrics=[tf.keras.metrics.MeanAbsoluteError()],
        )

        self.batch_size = params_dict["batch_size"]
        self.parton_data, self.reco_data = data_utils.load_jet_data(params_dict["data_path"])
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
            self.parton_data,
            self.reco_data,
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
        file_utils.save_network(self.save_dir, model_path="./FCNN.py")

    def save_params(self, params_dict):
        file_utils.save_params(self.save_dir, params_dict)