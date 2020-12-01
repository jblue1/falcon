import tensorflow as tf
from tensorflow import keras
import numpy as np
from tensorflow.keras.constraints import max_norm


class weight_clipping(tf.keras.constraints.Constraint):
    def __init__(self, ref_value):
        self.ref_value = ref_value



class cWGAN():
    def __init__(self, dataset, num_batches):
        # hyper parameters recommended by paper
        self.num_critic_iters = 5
        self.clip_value = 0.01
        self.critic_optimizer = tf.keras.optimizers.RMSprop(lr=5e-5)
        self.generator_optimizer = tf.keras.optimizers.RMSprop(lr=5e-5)
        self.batch_size = 64 

        self.generator = self.build_generator()
        self.critic = self.build_critic()
        self.dataset = dataset.batch(self.batch_size)
        self.dataset = self.dataset.shuffle(num_batches,
                reshuffle_each_iteration=True)
        


    def build_generator(self):
        noise = keras.Input(shape=(4,), name="noiseIn")
        pJet = keras.Input(shape=(4,), name="pjetIn")

        
        x = keras.layers.Dense(32, activation='relu', name="genx1")(pJet)
        x = keras.layers.Dense(32, activation='relu', name="genx2")(x)
        x = keras.layers.Dense(32, activation='relu', name="genx3")(x)

        y = keras.layers.Dense(32, activation='relu', name="geny1")(noise)
        y = keras.layers.Dense(32, activation='relu', name="geny2")(y)
        y = keras.layers.Dense(32, activation='relu', name="geny3")(y)

        concat = keras.layers.concatenate([x, y], name="concat")
        out = keras.layers.Dense(32, activation='relu', name="both1")(concat)
        out = keras.layers.Dense(32, activation='relu', name="both2")(concat)
        out = keras.layers.Dense(4, name="out")(out)
        return keras.Model([pJet, noise], out)
   

    def build_critic(self):
        pJet = keras.Input(shape=(4,))
        rJet = keras.Input(shape=(4,))

        x = keras.layers.Dense(32, activation='relu')(pJet)
        x = keras.layers.Dense(32, activation='relu')(x)
        x = keras.layers.Dense(32, activation='relu')(x)

        y = keras.layers.Dense(32, activation='relu')(rJet)
        y = keras.layers.Dense(32, activation='relu')(y)
        y = keras.layers.Dense(32, activation='relu')(y)

        concat = keras.layers.concatenate([x, y])
        out = keras.layers.Dense(32, activation='relu')(concat)
        out = keras.layers.Dense(32, activation='relu')(concat)
        out = keras.layers.Dense(1)(out)
        return keras.Model([pJet, rJet], out)


    def print_network(self):
        self.generator.summary()
        self.critic.summary()


    def critic_loss(self, real_output, fake_output):
        '''
        The negative of the estimate of the wasserstein distance (negative
        because we want to perform gradient ascent on the critic)
        '''
        loss = -(tf.math.reduce_mean(real_output) -
                tf.math.reduce_mean(fake_output))
        #print("Critic loss: {}".format(loss))
        return loss


    def generator_loss(self, fake_output):
        loss = -tf.math.reduce_mean(fake_output)
        #print("Gen Loss: {}".format(loss.numpy()))
        return loss

         
    def clip_critic_weights(self):
        for l in self.critic.layers:
            weights = l.get_weights()
            weights = [np.clip(w, -self.clip_value, self.clip_value) for w in weights]
            l.set_weights(weights)


    #@tf.function()
    def train_critic(self, pJets, rJets):
        shape = tf.shape(pJets)
        noise = tf.random.uniform(shape, 0, 1, tf.float32)
        with tf.GradientTape(persistent=True) as tape:
            generated_rJets = self.generator([pJets, noise],
                    training=False)
            real_output = self.critic([pJets, rJets], training=True)
            fake_output = self.critic([pJets, generated_rJets],
                    training=True)

            critic_loss = self.critic_loss(real_output, fake_output)

        
        critic_grads = tape.gradient(critic_loss,
                self.critic.trainable_variables)

        self.critic_optimizer.apply_gradients(zip(critic_grads,
            self.critic.trainable_variables))

        self.clip_critic_weights()


    #@tf.function()
    def train_generator(self, pJets):
        shape = tf.shape(pJets)
        noise = tf.random.uniform(shape, 0, 1, tf.float32)

        with tf.GradientTape() as tape:
            generated_rJets = self.generator([pJets, noise], training=True)
            fake_output = self.critic([pJets, generated_rJets], training=False)
            generator_loss = self.generator_loss(fake_output)

        generator_grads = tape.gradient(generator_loss,
                self.generator.trainable_variables)
        self.generator_optimizer.apply_gradients(zip(generator_grads,
            self.generator.trainable_variables))


    def train_step(self):
        data = self.dataset.take(self.num_critic_iters + 2)
        count = 0
        for batch in data:
            if count < self.num_critic_iters:
                self.train_critic(batch[0], batch[1])
            if count == self.num_critic_iters:
                self.train_generator(batch[0])
            if count == self.num_critic_iters + 1:
                pJets = batch[0]
                rJets = batch[1]
                #print("PJets")
                #print(pJets)
                #print("rJets")
                #print(rJets)
                shape = tf.shape(pJets)
                noise = tf.random.uniform(shape, 0, 1, tf.float32)

                generated_rJets = self.generator([pJets, noise], training=False)
                #print("Gen rJets")
                #print(generated_rJets)
                real_output = self.critic([pJets, rJets], training=False)
                fake_output = self.critic([pJets, generated_rJets], training=False)
                #print(real_output)
                #print(fake_output)
                wass_estimate = -self.critic_loss(real_output, fake_output)
            count += 1

        return wass_estimate.numpy()




def main():

    data = np.loadtxt('../../data/processed/matchedJets.txt', skiprows=2)
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

    data[:, 0] = (data[:, 0] - partonPtMin)/partonPtMax
    data[:, 1:3] = (data[:, 1:3] - partonMean)/partonStd
    data[:, 3] = (data[:, 3] - partonEMin)/partonEMax
    data[:, 4] = (data[:, 4] - pfPtMin)/pfPtMax
    data[:, 5:7] = (data[:, 5:7] - pfMean)/pfStd
    data[:, 7] = (data[:, 7] - pfEMin)/pfEMax

    np.random.shuffle(data)
    trainParton = data[:, :4]
    trainPf = data[:, 4:]
    train_dataset = tf.data.Dataset.from_tensor_slices((trainParton,
            trainPf))

    net = cWGAN(train_dataset, 100)
    net.clip_critic_weights()


if __name__ == "__main__":
    main()
