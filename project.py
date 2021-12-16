#Import library yg dibutuhkan
import os
import warnings
import numpy as np
import pandas as pd
import streamlit as st
from random import sample
from matplotlib import pyplot as plt
from matplotlib import colors as pltc

#Setting default
pd.set_option('display.float_format', lambda x: '%.2f' % x)
warnings.simplefilter(action='ignore', category=FutureWarning)

def getData():
    os.path.dirname(os.path.realpath('__file__'))
    #Baca data produksi minyak tiap negara
    dataprodminyak = pd.read_csv('produksi_minyak_mentah.csv')

    #Baca data kode negara lalu lakukan penyesuaian nama kolom
    datakodenegara = pd.read_json('kode_negara_lengkap.json')
    datakodenegara = datakodenegara.rename(columns={
                               'name': 'nama_negara',
                               'alpha-3': 'kode_negara',
                               'sub-region': 'subregion'})
    datakodenegara['nama_negara'] = datakodenegara['nama_negara'].str.replace(r'\(.+','')
    datakodenegara = datakodenegara[['nama_negara', 'kode_negara', 'region', 'subregion']]

    #Lakukan proses penggabungan data minyak dan kode negara
    data = dataprodminyak.merge(datakodenegara, how = 'inner', on = ['kode_negara'])
    
    #Urutkan berdasarkan nama_negara dan tahu6n produksi
    data = data.sort_values(by = ['nama_negara', 'tahun']).round(2).reset_index(drop = True)

    return data

def Konfigurasi(data):
    tahun  = sorted(data['tahun'].unique().tolist())
    negara = data['nama_negara'].unique().tolist()
    return tahun, negara

def ProduksiOptimal(data, namanegara):
    datanegara = data[data['nama_negara'] == namanegara]
    maxminprod = datanegara[(datanegara['produksi'].max() == datanegara['produksi']) | 
                            (datanegara['produksi'].min() == datanegara['produksi'])]\
                            .reset_index().sort_values(by = ['produksi'])\
                            .reset_index(drop = True)\
                            .drop_duplicates(subset=['produksi'], keep = 'last')

    return maxminprod

def plotProdMinyakTiapNegara(data, namanegara):
    #Inisialisasi data
    dataplot     = data[data['nama_negara'] == namanegara]
    tahunoptimum = ProduksiOptimal(dataplot, namanegara)['tahun'].tolist()
    prodoptimum  = ProduksiOptimal(dataplot, namanegara)['produksi'].tolist()

    #Inisialisasi absis dan ordinat
    x, y = dataplot['tahun'], dataplot['produksi']

    fig, ax = plt.subplots(1, figsize=(15,8))
    ax.plot(x, y, color = '#fe8a6c', linestyle='-', marker='o')  

    #Hapus frame line
    for spine in ['top','right']:
        ax.spines[spine].set_visible(False)
    
    #Atur batasan nilai y
    if((prodoptimum[0] > 0) & (prodoptimum[0] < 5000)):
        ax.set_ylim(bottom = -1000, top = 1.2*prodoptimum[1])
    elif(prodoptimum[0] == 0):
        ax.set_ylim(bottom = -100, top = 1.2*prodoptimum[1] + 100)
    else:
        ax.set_ylim(bottom = 0.7*prodoptimum[0], top = 1.2*prodoptimum[1])
    
    #Beri atribut pada grafik
    ax.set_title('Jumlah Produksi Minyak Negara {}\nPada Tahun {} - {}'
                  .format(namanegara, dataplot['tahun'].min(), dataplot['tahun'].max()),
                  color = '#f7374b',
                  fontsize = 15)
    
    ax.set_xticks(x)
    ax.set_xticklabels(x, rotation=55)
    ax.legend(['Produksi Minyak'])
    ax.set_xlabel("Tahun", labelpad = 10)
    ax.set_ylabel("Jumlah Produksi (barrel)", labelpad = 10)

    #Memberi anotasi pada data produksi maksimum dan minimum
    for x, y in zip(dataplot['tahun'], dataplot['produksi']):
        if ((x in tahunoptimum) & (y > 0)):
            label = "{:.2f}".format(y)
            if(x == tahunoptimum[0]):
                ax.annotate(label + "\n(min)", (x,y), textcoords = "offset points", xytext = (0, -30), ha='center')
            elif(x == tahunoptimum[1]):
                ax.annotate(label + "\n(max)", (x,y), textcoords = "offset points", xytext = (0, 10), ha='center')

    return fig, ax

def UrutanProduksiPerTahun(data, year, limit):
    datasorting = data[data['tahun'] == year]
    datasorting = datasorting.sort_values(by = ['produksi'], ascending = False)\
                            .reset_index().iloc[0 : limit]
    return datasorting

def plotUrutanProduksiPerTahun(data, year):    
    #Plot grafik
    colour = []
    for i in range(len(data['produksi'])):
        if(len(data['produksi']) >= 5):
            if(i < 3):
                colour.append('#77dd77')
            else:
                colour.append('#b0ffad')
        else:
            if(i == 0):
                colour.append('#77dd77')
            else:
                colour.append('#b0ffad')

    #Plot grafik
    x, y = data['nama_negara'], data['produksi']
    fig, ax = plt.subplots(figsize=(10,8))
    ax.barh(x, y, color = colour)

    #Invert sumbu y (max -> min)
    ax.invert_yaxis()

    ax.set_title("Urutan Negara Penghasil Minyak Terbanyak\npada Tahun {}".format(year), fontsize = 15, pad = 50)
    #ax.set_xlabel("Jumlah Produksi Minyak (barel)")
    #ax.set_ylabel("Negara")
    ax.text(1, 0.4, year, transform = ax.transAxes, color='#777777', size=46, ha='right', weight=800)
    ax.grid(which='major', axis='x', linestyle='--', alpha = 0.7)
    ax.xaxis.set_ticks_position('top')

    #Menghilangkan frame pada chart
    for i in ['bottom','left','right']:
      ax.spines[i].set_visible(False)

    #Memberi anotasi pada grafik
    j = 0
    for i in plt.gca().patches:
        ax.text(i.get_width()+.5, i.get_y()+.4, str(data['produksi'][j]), fontsize = 8, color='dimgrey')
        j = j + 1

    return fig, ax 

def UrutanProduksiKumulatif(data, limit):
    datasorting = data.groupby(['countryname'])['oilproduction'].sum().reset_index()
    datasorting = datasorting.sort_values(by = ['oilproduction'], ascending = False)\
                            .reset_index(drop = True).iloc[0:limit]
    return datasorting

class BubbleChart:
    def __init__(self, area, bubble_spacing=0):
        area = np.asarray(area)
        r = np.sqrt(area / np.pi)

        self.bubble_spacing = bubble_spacing
        self.bubbles = np.ones((len(area), 4))
        self.bubbles[:, 2] = r
        self.bubbles[:, 3] = area
        self.maxstep = 2 * self.bubbles[:, 2].max() + self.bubble_spacing
        self.step_dist = self.maxstep / 2

        length = np.ceil(np.sqrt(len(self.bubbles)))
        grid = np.arange(length) * self.maxstep
        gx, gy = np.meshgrid(grid, grid)
        self.bubbles[:, 0] = gx.flatten()[:len(self.bubbles)]
        self.bubbles[:, 1] = gy.flatten()[:len(self.bubbles)]

        self.com = self.center_of_mass()

    def center_of_mass(self):
        return np.average(self.bubbles[:, :2], axis=0, weights=self.bubbles[:, 3])

    def center_distance(self, bubble, bubbles):
        return np.hypot(bubble[0] - bubbles[:, 0], bubble[1] - bubbles[:, 1])

    def outline_distance(self, bubble, bubbles):
        center_distance = self.center_distance(bubble, bubbles)
        return center_distance - bubble[2] - \
            bubbles[:, 2] - self.bubble_spacing

    def check_collisions(self, bubble, bubbles):
        distance = self.outline_distance(bubble, bubbles)
        return len(distance[distance < 0])

    def collides_with(self, bubble, bubbles):
        distance = self.outline_distance(bubble, bubbles)
        idx_min = np.argmin(distance)
        return idx_min if type(idx_min) == np.ndarray else [idx_min]

    def collapse(self, n_iterations=50):
        for _i in range(n_iterations):
            moves = 0
            for i in range(len(self.bubbles)):
                rest_bub = np.delete(self.bubbles, i, 0)
                dir_vec = self.com - self.bubbles[i, :2]
                dir_vec = dir_vec / np.sqrt(dir_vec.dot(dir_vec))

                new_point = self.bubbles[i, :2] + dir_vec * self.step_dist
                new_bubble = np.append(new_point, self.bubbles[i, 2:4])

                if not self.check_collisions(new_bubble, rest_bub):
                    self.bubbles[i, :] = new_bubble
                    self.com = self.center_of_mass()
                    moves += 1
                else:
                    for colliding in self.collides_with(new_bubble, rest_bub):
                        dir_vec = rest_bub[colliding, :2] - self.bubbles[i, :2]
                        dir_vec = dir_vec / np.sqrt(dir_vec.dot(dir_vec))
                        orth = np.array([dir_vec[1], -dir_vec[0]])
                        new_point1 = (self.bubbles[i, :2] + orth *
                                      self.step_dist)
                        new_point2 = (self.bubbles[i, :2] - orth *
                                      self.step_dist)
                        dist1 = self.center_distance(
                            self.com, np.array([new_point1]))
                        dist2 = self.center_distance(
                            self.com, np.array([new_point2]))
                        new_point = new_point1 if dist1 < dist2 else new_point2
                        new_bubble = np.append(new_point, self.bubbles[i, 2:4])
                        if not self.check_collisions(new_bubble, rest_bub):
                            self.bubbles[i, :] = new_bubble
                            self.com = self.center_of_mass()

            if moves / len(self.bubbles) < 0.1:
                self.step_dist = self.step_dist / 2

    def plot(self, ax, labels, colors):
        for i in range(len(self.bubbles)):
            circ = plt.Circle(self.bubbles[i, :2], self.bubbles[i, 2], color=colors[i])
            ax.add_patch(circ)
            ax.text(*self.bubbles[i, :2], labels[i], horizontalalignment='center', verticalalignment='center', fontsize = 12)

def randomcolor(n):
    all_colors = [k for k, v in pltc.cnames.items()]
    for i in ['black', 'white', 'cyan', 'aqua']:
      all_colors.remove(i)
    colors = sample(all_colors, n)
    return colors

def plotBubbleChart(data, limit = 0, years = (0,0)):
    data['country_name'] = (data.index+1).astype('str') + '. ' + data['nama_negara'] + '\n' + data['produksi'].round(2).astype('str') 
    
    data['country_name'] = data['country_name'].str.replace('and ','\nand ')
    data['country_name'] = data['country_name'].str.replace('of ','of\n')

    data = data.sort_values(by = ['nama_negara']).reset_index()
    colors = randomcolor(len(data['nama_negara']))
    
    bubble_chart = BubbleChart(area = data['produksi'], bubble_spacing = 0.1)
    bubble_chart.collapse()

    fig, ax = plt.subplots(subplot_kw = dict(aspect = "equal"))
    fig.set_size_inches(15, 10)
    bubble_chart.plot(ax, data['country_name'], colors)
    ax.axis("off")
    ax.relim()
    ax.autoscale_view()
    ax.set_title('{} Negara Penghasil Minyak Terbanyak\n(Kumulatif pada Tahun {} - {})'.format(limit, years[0], years[1]), fontsize = 25)

    return fig, ax

def dataSummary(data):
    dataproduksinonzero = data[data['produksi'] > 0]
    dataproduksizero    = data[data['produksi'] == 0].sort_values(by = ['nama_negara', 'tahun']).reset_index(drop = True)
    maxdataproduksi     = dataproduksinonzero[dataproduksinonzero.groupby(['tahun'])['produksi'].transform(max) == dataproduksinonzero['produksi']]\
                                .sort_values(by = ['tahun']).reset_index(drop = True)
    mindataproduksi     = dataproduksinonzero[dataproduksinonzero.groupby(['tahun'])['produksi'].transform(min) == dataproduksinonzero['produksi']]\
                                .sort_values(by = ['tahun']).reset_index(drop = True)
    return maxdataproduksi, mindataproduksi, dataproduksizero

#Program Utama
data          = getData()
tahun, negara = Konfigurasi(data)

st.set_page_config(layout = "wide") 
judul1 = "Statistik Produksi Minyak di Dunia pada Tahun {} - {}".format(min(tahun), max(tahun)) 
st.title(judul1)
st.markdown("*Anggie Fiorella*")

#image = Image.open('tj_logo.png')
#st.sidebar.image(image)

st.sidebar.title("Pengaturan")
left_col, right_col = st.columns((1,1))

st.sidebar.subheader("Pengaturan konfigurasi tampilan")
T = st.sidebar.selectbox("Pilih Tahun", tahun)
N = st.sidebar.selectbox("Pilih Negara", negara)
B = st.sidebar.number_input("Jumlah negara yang ditampilkan pada urutan", min_value = 1, max_value = len(negara), value = 10)

left_col.subheader("Tabel Representasi Data Produksi Minyak Dunia")
left_col.dataframe(data)

#**** SOAL 1 ****#
fig1, plot1 = plotProdMinyakTiapNegara(data, N)
judul1 = 'Jumlah Produksi Minyak Negara {} Pada Tahun {} - {}'.format(N, min(tahun), max(tahun))
right_col.subheader(judul1)
right_col.pyplot(fig1)

#**** SOAL 2 ****#
dataterurut   = UrutanProduksiPerTahun(data, year = T, limit = B)
fig2, plot2   = plotUrutanProduksiPerTahun(dataterurut, year = T)

judul2 = "Urutan Negara Penghasil Minyak Terbanyak pada Tahun {}".format(T)
left_col.subheader(judul2)
left_col.pyplot(fig2)

#**** SOAL 3 ****#
datakumulatif = UrutanProduksiKumulatif(data, limit = B)
fig3, plot3   = plotBubbleChart(datakumulatif, limit = B, years = (min(tahun), max(tahun)))
judul3 = '{} Negara Penghasil Minyak Terbanyak (Kumulatif pada Tahun {} - {})'.format(B, min(tahun), max(tahun))
right_col.subheader(judul3)
right_col.pyplot(fig3)

maxdataproduksi, mindataproduksi, data0minyak = dataSummary(data)
left_col.subheader('Jumlah Produksi Minyak Terbesar di Dunia dari Tahun ke Tahun')
left_col.dataframe(maxdataproduksi)
right_col.subheader('Jumlah Produksi Minyak Paling Rendah di Dunia dari Tahun ke Tahun')
right_col.dataframe(mindataproduksi)

left_col.subheader('Negara yang Tidak Menghasilkan Minyak Pada Tahun {} - {}'.format(min(tahun), max(tahun)))
left_col.dataframe(data0minyak)