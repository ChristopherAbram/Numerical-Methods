# Metody-Obliczeniowe

W tym repo znajdują się przykładowe implementacje podstawowych metod obliczeniowych. Poniżej treści zadań:

1. Napisz program w języku „C/C++”, umożliwiający „doświadczalne” wyznaczenie liczby bitów
mantysy oraz tzw. epsylona maszynowego, dla zmiennych typu float i double, tj. najmniejszej liczby
e takiej, że fl(e + 1) > 1. Jaki jest związek e z precyzją arytmetyki?

2. Zaimplementuj w języku „C/C++” algorytm obliczający przybliżone wartości funkcji f(x) = exp(x) dla x  [-30, 30], 
poprzez sumowanie N wyrazów rozwinięcia tej funkcji w szereg Taylora wokół x = 0. 
Zbadaj jak zmieniają się błędy względne przybliżenia funkcji w tym algorytmie, przy wzrastającej liczbie N  [1, 1000]. 
Wyjaśnij przyczyny obserwowanych błędów i ich zmian ze wzrostem N. 
Następnie dokonaj takiej modyfikacji algorytmu, i wybierz takie N, aby uzyskać dokładność maszynową dla dowolnego x  [-30, 30]. 
W obliczeniach zastosuj zmienne podwójnej precyzji.

3. Napisz program w języku „C/C++”, realizujący metody:
  (a) Picarda
  (b) bisekcji
  (c) Newtona
  (d) siecznych
rozwiązywania pojedynczych algebraicznych równań nieliniowych. Zastosuj program do
przykładów z zadania 1. Zastosuj trzy niezależne kryteria zakończenia iteracji. Zadbaj o to, aby
wyprowadzać na konsolę wyniki pośrednie obliczeń dla każdej iteracji, tak aby możliwe było
obserwowanie zbieżności kolejnych przybliżeń pierwiastków i porównanie liczby iteracji
niezbędnych do uzyskania rozwiązania o zadanej dokładności przez każdą z metod. 

  1) cos^2(x/4) - x = 0				-> f1(x) = cos^2(x/4) - x
  2) exp(x) - exp(-x) + x - 1 = 0		-> f2(x) = exp(x) - exp(-x) + x - 1

4. Napisz program w języku „C/C++”, realizujący metodę Newtona rozwiązywania układu trzech 
algebraicznych równań nieliniowych, i zastosuj ten program do przykładu z zadania 1. 
Przyjmij takie przybliżenie początkowe, aby uzyskać zbieżność metody. Zastosuj trzy niezależne kryteria 
zakończenia iteracji. Zadbaj o to, aby wyprowadzać na konsolę wyniki pośrednie obliczeń dla każdej iteracji, 
tak aby możliwe było obserwowanie zbieżności kolejnych przybliżeń pierwiastków i porównanie liczby iteracji 
niezbędnych do uzyskania rozwiązania o zadanej dokładności. 
Oblicz jak zmienia się residuum układu w trakcie iteracji.

  (1)   xy - 2 = 0
  (2)   y/2 - sin(pi/4 - z) = 0
  (3)   x^2 + y^2 + z^2 - 4 = 0
  
5. Dana jest macierz A:

	1   -20   30  -4
	2   -40  -6    50
	9   -180  11  -12
   -16   15  -140  13

oraz wektor b:

	 35
	 104
	-366
	-354

Napisz program w języku „C/C++”, realizujący dekompozycję LU macierzy A, przy zastosowaniu
eliminacji Gaussa z częściowym wyborem elementu podstawowego, a następnie rozwiązujący układ
równań Ax = b. Uwaga: należy zrealizować wariant dekompozycji omawiany na wykładzie.

6. Napisz program w języku „C/C++”, realizujący algorytm Thomasa dla macierzy trój-diagonalnej o
dowolnych rozmiarach N x N, a następnie zastosuj ten program do rozwiązania układu równań
Ax = b, w którym

Macierz A:
	10   1/2
	1/3  20   1/4
	     1/5  30   1/6
	          1/7  30   1/8
	               1/9  20    1/10
					            1/11  10
Macierz b:
	3
	165/4
	917/30
	851/28
	3637/90
	332/11

Program należy zrealizować w postaci dwóch oddzielnych funkcji: jednej, która korzysta wyłącznie
z danych dotyczących macierzy A, oraz drugiej, która korzysta z wektora b oraz wyników działania
pierwszej funkcji dla macierzy A.

7. Napisz program w języku „C/C++”, rozwiązujący układ czterech równań liniowych metodami
iteracyjnymi: 
	(a) Jacobiego, 
	(b) Gaussa-Seidela, 
	(c) SOR z parametrem w = 1/2, 
a następnie zastosuj ten program do rozwiązania układu równań liniowych Ax = b, gdzie:

macierz A:

	100 -1  2   -3
	1   200 -4  5
	-2  4   300 -6
	3   -5  6   400

macierz b:

	116
	-226
	912
	-1174

Przyjmij przybliżenie początkowe x0:

	2
	2
	2
	2

Zastosuj trzy niezależne kryteria zakończenia iteracji. Zadbaj o to, aby wyprowadzać na konsolę
wyniki pośrednie obliczeń dla każdej iteracji, tak aby możliwe było obserwowanie zbieżności
kolejnych przybliżeń pierwiastków i porównanie liczby iteracji niezbędnych do uzyskania
rozwiązania o zadanej dokładności bezwzględnej. Oblicz jak zmienia się residuum układu w trakcie
kolejnych iteracji.

8. Napisz program w języku „C/C++”, obliczający przybliżone wartości pierwszych pochodnych
funkcji f(x) = cos(x) w punktach końcowych i środkowym przedziału [0, pi/2] zmiennej x. Zastosuj
wszystkie omawiane na wykładzie i na ćwiczeniach przybliżenia różnicowe dwupunktowe i
trzypunktowe (jednostronne bądź centralne, w zależności od położenia punktu w przedziale) na sieci
jednorodnej o kroku h. 

Wykonaj (na jednym rysunku) wykresy przedstawiające zależności błędów
bezwzględnych przybliżeń różnicowych od kroku sieci, posługując się skalą logarytmiczną (tzn.
wykresy zależności log10|błędu| od log10 h ). Na podstawie wykresów wyznacz doświadczalnie rzędy
dokładności przybliżeń różnicowych. Sprawdź, czy tak wyznaczone rzędy dokładności pokrywają
się z rzędami teoretycznymi i wyjaśnij ewentualne rozbieżności. Ponadto zidentyfikuj wartości
kroku sieci poniżej których pojawia się wpływ błędów maszynowych. Obliczenia powtórz dla
dwóch typów zmiennych rzeczywistych (float, i double) i porównaj wyniki.

Uwaga: najwygodniej jest zastosować wzorzec funkcji (function template) z typem zmiennych jako
parametrem wzorca.

9. Napisz program w języku „C/C++”, rozwiązujący równanie różniczkowe zwyczajne drugiego rzędu:
U''(x) + U(x) + 2sin(x) = 0, określone na przedziale 0 <= x <= pi/2,
z warunkami brzegowymi U(0) = 0, U(pi/2) = 0. Zastosuj typ double oraz trzypunktową
dyskretyzację konwencjonalną oraz dyskretyzację Numerowa na sieci jednorodnej. Do rozwiązania
układu liniowych równań algebraicznych zastosuj algorytm Thomasa (patrz zajęcia nr 6). Wykonaj
rusynek przedstawiający porównanie uzyskanych wyników numerycznych z rozwiązaniem
analitycznym U(x) = xcos(x) . Pokaż, że rząd dokładności
rozwiązań numerycznych jest zgodny z przewidywaniami teoretycznymi wynikającymi z zadania 2.
W tym celu wykonaj (na jednym rysunku) wykresy przedstawiające zależności maksymalnego błędu
bezwzględnego rozwiązań od kroku sieci h, posługując się skalą logarytmiczną (tzn. wykresy
zależności log10|błędu| od log10 h ). Na podstawie wykresów wyznacz doświadczalnie rzędy
dokładności rozwiązań uzyskanych za pomocą różnych metod, i porównaj je z rzędami
teoretycznymi. Ponadto zidentyfikuj wartości kroku sieci poniżej których pojawia się wpływ błędów
maszynowych.

10. Napisz program w języku „C/C++”, rozwiązujący równanie różniczkowe zwyczajne pierwszego rzędu:
	dy(t)/dt + (10t^2 + 20)/(t^2 + 1) * [y(t) - 1] = 0, określone dla zmiennej t >= 0,
z warunkiem początkowym y(0) = 0, za pomocą metod:

	(a) bezpośredniej Eulera,
	(b) pośredniej Eulera,
	(c) metody trapezów.

Dla metod (b) i (c) wykonaj oddzielne rysunki przedstawiające po dwa wykresy: 
wykres przykładowego rozwiązania numerycznego oraz (dla porównania) wykres rozwiązania
analitycznego: y(t) = 1 - exp{-10[t + arctg(t)]}. Oba wykresy winny przedstawiać zależność y od
zmiennej niezależnej t. Rozwiązania analityczne zaznacz linią ciągłą, a numeryczne punktami. W
przypadku metody (a) wykonaj dwa takie rysunki: jeden uzyskany w warunkach numerycznej
stabilności metody, a drugi w warunkach numerycznej niestabilności. Wyjaśnij różnice pomiędzy
uzyskanymi wykresami.

Pokaż, że rząd dokładności uzyskanych stabilnych rozwiązań numerycznych jest zgodny z
przewidywaniami teoretycznymi. W tym celu wykonaj (na jednym rysunku) wykresy
przedstawiające zależności maksymalnych błędów bezwzględnych rozwiązań uzyskanych trzema
metodami, od kroku sieci czasowej (delta)t, posługując się skalą logarytmiczną (tzn. wykresy zależności
log10|błędu| od log10(delta)t ). Na podstawie wykresów wyznacz doświadczalnie rzędy dokładności
rozwiązań uzyskanych za pomocą różnych metod i porównaj je z rzędami teoretycznymi. O ile to
możliwe, zidentyfikuj też wartości kroku sieci poniżej których pojawia się wpływ błędów
maszynowych.

11. Napisz program w języku „C/C++”, demonstrujący zjawisko Rungego w interpolacji wielomianowej Lagrange’a, 
na przykładzie interpolacji funkcji f(x) = x/(1 + 25|x|^3), określonej w przedziale [-1, 1]. 
Porównaj wyniki interpolacji na węzłach równoodległych z wynikami interpolacji na węzłach Czebyszewa.
