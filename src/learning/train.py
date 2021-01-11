import sys
import file_utils
import data_utils
import FCNN
import cWGAN
import tensorflow as tf
import os



class FCNNTrainer():
    def __init__(self, lr, epochs, data_path, batch_size):
        self.save_dir = file_utils.make_save_directory("FCNN")
        self.model = FCNN.make_model()
        self.model.compile(optimizer=tf.keras.optimizers.Adam(learning_rate=lr),
           loss=tf.keras.losses.MeanAbsoluteError(),
           metrics=[tf.keras.metrics.MeanAbsoluteError()])

        self.batch_size = batch_size
        self.parton_data, self.reco_data = data_utils.load_jet_data(data_path)
        self.checkpoint_path = os.path.join(self.save_dir, 'training/cp.cpkt')
        self.epochs = epochs


    def train(self):
        cp_callback = tf.keras.callbacks.ModelCheckpoint(filepath=self.checkpoint_path,
            save_weights_only=True,
            verbose=1,
            save_best_only=True)

        history = self.model.fit(self.parton_data,
            self.reco_data,
            batch_size=self.batch_size,
            epochs=self.epochs,
            validation_split=0.2,
            callbacks=[cp_callback])

        return history


    def save_losses(self, history):
        training_loss = history.history['loss']
        val_loss = history.history['val_loss']
        loss_dict = {'Training Loss' : training_loss, 'Validation Loss' : val_loss}
        file_utils.save_losses(self.save_dir, loss_dict)


    def save_model(self):
        file_utils.save_network(self.save_dir, model_path="./FCNN.py")


def train_fcnn(params):
    print("Training fcnn")
    trainer = FCNNTrainer(lr=1e-4, epochs=1, data_path='../../data/processed/newPartonMatchedJets.txt', batch_size=64)
    history = trainer.train()
    trainer.save_losses(history)
    trainer.save_model()


def train_cGAN(params):
    print("Training cGAN")


def train_cWGAN(params):
    print("Training cWGAN")


def train(model, params):
    if model == "FCNN":
        train_fcnn(params)

    elif model == "cGAN":
        train_cGAN(params)

    elif model == "cWGAN":
        train_cWGAN(params)

    else:
        usage("Invalid model chosen")


def usage(message):
    print(message)
    print("usage: python train.py model param_file")
    print("    model - Model to train. One of FCNN, cGAN, cWGAN")
    print("    param_file - path to json file with model parameters")
    sys.exit(1)


def main():
    if len(sys.argv) != 3:
        usage("Incorrect number of arguments given")
    model = sys.argv[1]
    params = sys.argv[2]
    train(model, params)


if __name__ == "__main__":
    main()