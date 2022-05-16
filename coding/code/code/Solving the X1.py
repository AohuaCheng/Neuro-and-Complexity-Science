import tensorflow as tf
from tensorflow import keras
import numpy as np

class MLP(kf.keras.model):
    def __init__(self):
        super().__init__()
        self.flatten = tf.keras.layers.Flatten()
        self.dense1 = tf.keras.layers.Dense(units=1000, activation=tf.nn.relu)
