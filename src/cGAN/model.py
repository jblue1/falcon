import tensorflow as tf
from tensorflow import keras


def make_generator():
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
   

def make_discriminator():
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
    out = keras.layers.Dense(1, activation='sigmoid')(out)
    return keras.Model([pJet, rJet], out)

def main():
    gen = make_simple_gen()
    gen.summary()

if __name__ == "__main__":
    main()
