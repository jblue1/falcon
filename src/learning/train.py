import sys
import file_utils
import data_utils
import FCNN
import cWGAN
import tensorflow as tf
import os
import time
import numpy as np


class FCNNTrainer:
    def __init__(self, lr, epochs, data_path, batch_size):
        self.save_dir = file_utils.make_save_directory("FCNN")
        self.model = FCNN.make_model()
        self.model.compile(
            optimizer=tf.keras.optimizers.Adam(learning_rate=lr),
            loss=tf.keras.losses.MeanAbsoluteError(),
            metrics=[tf.keras.metrics.MeanAbsoluteError()],
        )

        self.batch_size = batch_size
        self.parton_data, self.reco_data = data_utils.load_jet_data(data_path)
        self.checkpoint_path = os.path.join(self.save_dir, "training/cp.cpkt")
        self.epochs = epochs

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


class cWGANTrainer:
    def __init__(self, params_dict):
        self.save_dir = file_utils.make_save_directory("cWGAN")
        self.num_critic_iters = params_dict["num_critic_iters"]
        self.batch_size = params_dict["batch_size"]
        clip_value = params_dict["clip_value"]
        noise_dims = params_dict["noise_dims"]
        self.model = cWGAN.cWGAN_mnist(clip_value, noise_dims)

        self.data = data_utils.load_mnist_data()
        self.epochs = params_dict["epochs"]
        self.num_training_examples = len(self.data[0])
        self.weight_saving_interval = params_dict["weight_saving_interval"]

    def sample_batch_of_data(self):
        indices = np.random.choice(
            np.arange(self.num_training_examples), self.batch_size, replace=False
        )
        labels = self.data[1][indices]
        images = self.data[2][indices]
        return labels, images

    def train(self):
        """
        Training loop for the cWGAN. An "epoch" is considered to be when the generator
        has seen the same number of examples as are in the data set (note that since the
        batches are randomly sampled, it won't actually get trained on all the data
        each epoch).
        """
        critic_losses = []
        generator_losses = []
        wass_estimates = []
        batches_per_epoch = self.num_training_examples // self.batch_size
        for epoch in range(self.epochs):
            start = time.time()
            for batch_number in range(batches_per_epoch):
                # train critic for num_critic_iters
                for critic_iter in range(self.num_critic_iters):
                    labels, images = self.sample_batch_of_data()
                    critic_loss = self.model.train_critic(labels, images)
                    critic_losses.append(critic_loss)
                    self.model.clip_critic_weights()
                # train generator
                labels, images = self.sample_batch_of_data()
                concat_real = data_utils.concatenate_images_labels(images, labels)
                generator_loss = self.model.train_generator(labels)
                predicted_images = self.model.make_generator_predictions(labels)
                concat_fake = data_utils.concatenate_images_labels(
                    predicted_images, labels
                )
                generator_losses.append(generator_loss)
                real_output = self.model.critic(concat_real, training=False)
                fake_output = self.model.critic(concat_fake, training=False)
                wass_estimate = -self.model.critic_loss(real_output, fake_output)
                wass_estimates.append(wass_estimate)
                iteration = epoch * batches_per_epoch + batch_number
                if iteration % self.weight_saving_interval == 0:
                    self.save_weights(iteration)
                print(
                    "Iteration: {}  Wasserstein Estimate: {}".format(
                        iteration, wass_estimate
                    )
                )
            print("Time for epoch {}: {:1f}s".format(epoch, time.time() - start))
        return critic_losses, generator_losses, wass_estimates

    def save_weights(self, iteration):
        checkpoint_dir = self.save_dir + "/training_checkpoints"
        gen_filename = os.path.join(checkpoint_dir, "gen_" + str(iteration))
        critic_filename = os.path.join(checkpoint_dir, "critic_" + str(iteration))
        print("Saving generator weights at {}".format(gen_filename))
        self.model.generator.save_weights(gen_filename)
        self.model.critic.save_weights(critic_filename)

    def save_losses(self, critic_losses, generator_losses, wass_estimates):
        critic_loss_dict = {
            "Critic Loss": critic_losses,
        }
        file_utils.save_losses(self.save_dir, critic_loss_dict, "critic_")
        generator_loss_dict = {
            "Generator Loss": generator_losses,
            "Wasserstein Estimates": wass_estimates,
        }
        file_utils.save_losses(self.save_dir, generator_loss_dict, "generator_")

    def save_model(self):
        file_utils.save_network(self.save_dir, model_path="./cWGAN.py")

    def save_params(self, params_dict):
        file_utils.save_params(self.save_dir, params_dict)


def train_fcnn(params):
    print("Training fcnn")
    trainer = FCNNTrainer(
        lr=1e-4,
        epochs=1,
        data_path="../../data/processed/newPartonMatchedJets.txt",
        batch_size=64,
    )
    history = trainer.train()
    trainer.save_losses(history)
    trainer.save_model()


def train_cGAN(params):
    print("Training cGAN")


def train_cWGAN(params):
    print("Training cWGAN")
    params_dict = file_utils.get_cWGAN_hyperparams(params)
    trainer = cWGANTrainer(params_dict)
    critic_losses, generator_losses, wass_estimates = trainer.train()
    trainer.save_losses(critic_losses, generator_losses, wass_estimates)
    trainer.save_model()
    trainer.save_params(params_dict)


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