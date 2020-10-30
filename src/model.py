import tensorflow as tf
from tensorflow import keras

def make_model():
    model = keras.Sequential()
    model.add(keras.Input(shape=(4,)))
    model.add(keras.layers.Dense(128, activation='relu'))
    model.add(keras.layers.Dense(128, activation='relu'))
    model.add(keras.layers.Dense(128, activation='relu'))
    model.add(keras.layers.Dense(128, activation='relu'))
    model.add(keras.layers.Dense(128, activation='relu'))
    model.add(keras.layers.Dense(4))
    return model
