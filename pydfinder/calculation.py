import numpy as np
from sklearn import metrics
from sklearn.cluster import KMeans
import scipy.stats as stats


def partition(input_list):
    sorted_list = sorted(input_list)
    max_index = np.diff(sorted_list).argmax()
    list1 = sorted_list[:max_index+1]
    list2 = sorted_list[max_index+1:]
    return list1, list2


def silhouette_score(list1, list2):
    X = np.array(list1 + list2)
    y = np.array([0]*len(list1) + [1]*len(list2))
    return metrics.silhouette_score(X.reshape(-1, 1), y)

def perform_f_test(input_list):
    l1, l2 = partition(input_list)
    print l1, l2
    return stats.f_oneway(l1, l2)

def main2():
    print perform_f_test(
        [532, 570, 470, 486, 563, 571, 383, 463, 2462, 551, 502, 484, 448, 586, 603, 457, 568, 512, 402, 495]
    )
    print perform_f_test(
        [387, 555, 550, 517, 506, 449, 490, 511, 487, 392]
    )
    print perform_f_test(
        [474, 493, 390, 457, 481, 445, 484, 583, 460, 506, 430]
    )


def main():
    input_list = [532, 570, 470, 486, 563, 571, 383, 463, 2462, 551, 502, 484, 448, 586, 603, 457, 568, 512, 402, 495]
    X = np.array(input_list).reshape(-1, 1)
    km = KMeans(2)
    y_km = km.fit_predict(X)
    print metrics.silhouette_score(X, y_km)
    l1, l2 = partition(input_list)
    print silhouette_score(l1, l2)
    input_list = [387, 555, 550, 517, 506, 449, 490, 511, 487, 392]
    X = np.array(input_list).reshape(-1, 1)
    km = KMeans(2)
    y_km = km.fit_predict(X)
    print y_km
    print metrics.silhouette_score(X, y_km)
    l1, l2 = partition(input_list)
    print l1, l2
    print silhouette_score(l1, l2)
    input_list = [474, 493, 390, 457, 481, 445, 484, 583, 460, 506, 430]
    X = np.array(input_list).reshape(-1, 1)
    km = KMeans(2)
    y_km = km.fit_predict(X)
    print y_km
    print metrics.silhouette_score(X, y_km)
    l1, l2 = partition(input_list)
    print l1, l2
    print silhouette_score(l1, l2)

if __name__ == '__main__':
    # main()
    main2()
