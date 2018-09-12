#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
\author Pawel Czajka

do creatora dajemy sobie nazwe pliku (sciezke) oraz slownik typu {'czterowektor': (4,'f'),'intowa_wlasnosc': (1,'i'),
...}
to znaczy nazwe, ile to jest liczb, jakiego typu. Obsluguje na razie jedynie 'f' oraz 'i' to 
jest float oraz int
nie mozna uzyc jako klucz 'label', bo to jest wykorzystywana nazwa.


funkcja wpisz bierze jako argument generator ktory zwraca jak na niego podziałać next() cos typu 
({"czterowektor":[1.,2.,3.,4.], ...},1) gdzie 1 jest labelem, label jest intem.
. oczywiście klucze slownika zgadzają się 
z kluczami ze slownika ktorego uzylismy do kreatora.

dataset wczytany metoda wczytaj dataset jest juz w postaci wygodnej dla mnie to znaczy 
slownik feature, label
"""
import tensorflow as tf
import numpy as np
import matplotlib.pyplot as plt
from sklearn.utils import shuffle

class Io_tf_binary:
    def __init__(self,nazwa_pliku,slownik):
        self.plik=nazwa_pliku
        self.typy=slownik
        #self.cos=Io_tf_binary.wrap_int64([5])
    

    """
    To moze sobie czytac ktos kto chce zmieniac wnetrznosci tej klasy
    Nie polecam
    
    Teraz ta funkcja wpisz jest ważna. ona bierze po jednym przypadku testowym, 
    ( to jest ta petla for i in range()) i go zapisuje. trzeba zwracac 
    uwage na to jakiego typu sa zapisywane rzeczy. mozna oczywiscie zrobic slownik
    data dluzszym, jesli to w jakis sposob ulatwi nam myslenie o naszych danych. 
    Bo te nasze dane to bedzie slownik list, w ktorych to listach rzeczy maja 
    juz taki sam typ, a klucze to beda jakies opisowe nazwy.
    np 

    data={
    'czteroped_lewej_nogi_czy_cos': wrap_float64(cztero), # gdzie cztero to jest tensor floatow o shape (4,)
    # reszta rzeczy

    }

    Jak byscie chcieli jako wartosci miec stringi to musicie pomyslec jak zrobic wrapy dla stringow. oczywiscie
    nie znajdziecie zadnej dokumentacji.

    UWAGA 
    w tym slowniku data musi byc to co klasyfikujemy oznaczone przy pomocy 'label' bo inaczej sie  wywali program.
    """


    def wpisz(self,generator):
        """tworzy ten nasz dataset w pliku out_path """
        
        def wrap_int64(value):
            """lista intow musi wlesc"""
            return tf.train.Feature(int64_list=tf.train.Int64List(value=value))
        def wrap_float64(value):
            """lista floatow musi wlesc"""
            return tf.train.Feature(float_list=tf.train.FloatList(value=value))
        
        #f,l=kolko_w_kolku() #mozna zmienic jak sie podoba
        def data_slownik(f,l):
            wyrzut={}
            for k in self.typy.keys():
                if self.typy[k][1]=='f':
                    wyrzut[k]=wrap_float64(np.array(f[k]).reshape((-1,)))
                else:
                    wyrzut[k]=wrap_int64(np.array(f[k]).reshape((-1,)))
            wyrzut['label']=wrap_int64([l])
            return wyrzut
                    
        with tf.python_io.TFRecordWriter(self.plik) as writer:
            for f,l in generator:
                data=data_slownik(f,l)
                # Wrap the data as TensorFlow Features.
                feature = tf.train.Features(feature=data)

                # Wrap again as a TensorFlow Example.
                example = tf.train.Example(features=feature)

                # Serialize the data.
                serialized = example.SerializeToString()

                # Write the serialized data to the TFRecords file.
                writer.write(serialized)

            
    def wczytaj_dataset(self):
        
        def zeslownikoj(x):
            keys=list(x.keys())
            f={}
            for k in keys:
                if not k=='label':
                    f[k]=x[k]
            return f,x['label']
        def features_generoj():
            wyrzut={}
            for k in self.typy.keys():
                if self.typy[k][1]=='f':
                    wyrzut[k]=tf.FixedLenFeature([self.typy[k][0]], tf.float32)
                else:
                    wyrzut[k]=tf.FixedLenFeature([self.typy[k][0]], tf.int64)
            wyrzut['label']=tf.FixedLenFeature([], tf.int64)
            return wyrzut

        def parse(serialized):
            
            # Define a dict with the data-names and types we expect to
            # find in the TFRecords file.
            # It is a bit awkward that this needs to be specified again,
            # because it could have been written in the header of the
            # TFRecords file instead.
            """
            features = \
                {
                    'dwuwektor': tf.FixedLenFeature([2], tf.float32),#z jakiegos powodu to jest float32, nie wiem czemu
                    'label': tf.FixedLenFeature([], tf.int64)
                }
            """
            features=features_generoj()
            print(features)

            # Parse the serialized data so we get a dict with our data.
            parsed_example = tf.parse_single_example(serialized=serialized,
                                                     features=features)


            return zeslownikoj(parsed_example)

        dataset = tf.data.TFRecordDataset(self.plik)
        dataset = dataset.map(parse)
        return dataset
    
    
    
    
