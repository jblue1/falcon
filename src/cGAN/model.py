import tensorflow as tf
from tensorflow import keras

def make_generator():
    noise = keras.Input(shape=(4,))
    pJet = keras.Input(shape=(4,))

    
    x = keras.layers.Dense(32, activation='relu')(pJet)
    x = keras.layers.Dense(32, activation='relu')(x)

    y = keras.layers.Dense(32, activation='relu')(noise)
    y = keras.layers.Dense(32, activation='relu')(y)

    concat = keras.layers.concatenate([x, y])
    out = keras.layers.Dense(32, activation='relu')(concat)
    out = keras.layers.Dense(4)(out)
    return keras.Model([pJet, noise], out)
   

def make_discriminator():
    pJet = keras.Input(shape=(4,))
    rJet = keras.Input(shape=(4,))

    x = keras.layers.Dense(32, activation='relu')(pJet)
    x = keras.layers.Dense(32, activation='relu')(x)

    y = keras.layers.Dense(32, activation='relu')(pJet)
    y = keras.layers.Dense(32, activation='relu')(y)

    concat = keras.layers.concatenate([x, y])
    out = keras.layers.Dense(32, activation='relu')(concat)
    out = keras.layers.Dense(1, activation='sigmoid')(out)
    return keras.Model([pJet, rJet], out)

