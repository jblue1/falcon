import sys
import file_utils
import data_utils
import FCNN
import cWGAN
import tensorflow as tf
import os
import time
import numpy as np


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


class cWGANTrainer():
    def __init__(self, data_path, num_critic_iters, batch_size, lr, epochs):
        self.save_dir = file_utils.make_save_directory("cWGAN")
        self.num_critic_iters = num_critic_iters
        self.batch_size = batch_size
        self.model = cWGAN.cWGAN_mnist(num_critic_iters, batch_size)
        self.data = data_utils.load_mnist_data()
        self.epochs = epochs
        self.num_training_examples = len(self.data[0])


    def train(self):
        batches_per_epoch = self.num_training_examples // self.batch_size 
        for epoch in range(self.epochs):
            start = time.time()

            for batch_number in range(5): #range(batches_per_epoch):

                # train critic for num_critic_iters
                for critic_iter in range(self.num_critic_iters):
                    indices = np.random.choice(np.arange(self.num_training_examples), self.batch_size, replace=False)
                    labels = self.data[1][indices]
                    images = self.data[2][indices]

                    critic_loss = self.model.train_critic(labels, images)
                    self.model.clip_critic_weights()
                    print("Critic loss {}".format(critic_loss))
                
                # train generator
                indices = np.random.choice(np.arange(self.num_training_examples), self.batch_size, replace=False)
                labels = self.data[1][indices]
                self.model.train_generator(labels)
                





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
    trainer = cWGANTrainer(data_path="../../data/processed/newPartonMatchedJets.txt", num_critic_iters=5,
    batch_size=64, lr=1e-5, epochs=1)
    trainer.train()


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