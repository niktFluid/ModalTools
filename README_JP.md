# ModalTools
圧縮性時間平均流れ場に関する線形Navier-Stokesオペレータの構築とモード解析用のコードです．
(Click [here](README.md) for English version.)

## 概要
このコードでは，線形化Navier-Stokes (LNS) オペレータの構築とその分析 (モード解析) を行うことができます．プロジェクトにはLNS方程式の右辺を評価するサブルーチンと，行列表現されたLNSオペレータを生成するコード (いわゆるzero-oneトリックを使用しています) ，および全体安定性解析やレゾルベント解析 ([参考文献](https://doi.org/10.2514/1.J056060)) 用のコードが含まれています．プロジェクトには，デモンストレーション用に[OpenFOAM](https://www.openfoam.com/)による[2次元円柱周り流れのシミュレーションデータ](CylinderFlow)が付属していますので，これらの解析を手軽に試すことができます．


## デモ
2次元円柱周り流れのシミュレーションデータを用いた，本コードの簡単なデモンストレーションを紹介します．流れ場の気流条件は<img src="https://latex.codecogs.com/gif.latex?Re=150">, <img src="https://latex.codecogs.com/gif.latex?M_\infty=0.2">としています. シミュレーションはOpenFOAMを用いて実行しました. 上記の条件でシミュレーションを実施した場合，円柱背後にはカルマン渦列が形成され，その振動周波数はストローハル数で<img src="https://latex.codecogs.com/gif.latex?St=0.178">になります. 時間平均流れ場を用いてLNSオペレータを構築し，安定性解析とレゾルベント解析を実施します．

![Figure1-2](https://user-images.githubusercontent.com/47338366/76111253-cc8fdc00-5f94-11ea-880f-b9a0953f5330.png)


### 安定性解析
流れ場の (全体) 安定性解析はLNSオペレータに関する固有値問題を解くことで実施できます．LNSオペレータのスペクトル (固有値の集合) は時間平均流れ場における微小擾乱の動的な成長特性に関する情報を，対応する固有ベクトルは流れ場の構造 (モード) の情報を持っています．固有値の実部はモードの成長率を表しており，虚部は周波数に対応しています．円柱周り流れの安定性解析の結果を下図に示します．スペクトル (下図左) を見ると，いくつかの固有値が不安定面 (<img src="https://latex.codecogs.com/gif.latex?{\rm&space;Re}(\lambda)&space;>&space;0">) に位置していることが分かります．この内，より高周波のモードに対応している固有値 (<img src="https://latex.codecogs.com/gif.latex?\lambda&space;=&space;7.25&space;\times&space;10^{-3}&space;&plus;&space;0.208i">) の虚部の値はストローハル数換算の周波数で<img src="https://latex.codecogs.com/gif.latex?St=0.165">に対応しており，シミュレーションで得られた渦振動の周波数に近い値を示しています．実際に，この固有値に対応する流れ場のモードは，下図右に示していますように，明らかに円柱背後の渦振動に対応しています．また，実軸上にある不安定モードに対応する固有値は，シミュレーションにおける時間平均操作が十分でないために現れていると考えられます (これらの固有値に対応するモードを可視化してみて下さい)．

![Figure2-2](https://user-images.githubusercontent.com/47338366/76111278-d580ad80-5f94-11ea-966a-903eed397f2a.png)


### レゾルベント解析
レゾルベント解析では，LNSオペレータで表現される線形システムの入力 (加振) ・出力 (応答) モードを調べることができます．LNSオペレータの特異値分解 (SVD: Singular Value Decomposition) により得られる左右特異行列に含まれる各特異ベクトル (モード) はそれぞれ線形システムの出力・入力空間の直交基底を構成しており，対応する特異値は入出力モード間のエネルギー増幅率 (ゲイン) を表しています．レゾルベント解析ではそれぞれの周波数に対して複数のモードが得られます．通常はSVDによって全てのモードを計算することは無く (膨大な時間とメモリを消費します)，特異値 (ゲイン) の大きい順に少数のモードに限って計算します．最大特異値に対応するモードは支配モードと呼ばれます．レゾルベント解析によって得られた，円柱周り流れに対するゲインの内，最も大きい3つのモードに対応するゲインの分布を下図左に示します．支配モードのゲインは渦振動の周波数付近でピークに達していますが，これはカルマン渦振動が流体力学的な不安定性に起因していることが原因です (入出力間のエネルギー比が極めて大きくなり得る).．下図右には渦振動の周波数における支配モードの可視化結果を示していますが，出力モードは明らかに渦振動のモードに対応していることが分かります．一方，入力モードは円柱近傍における何らかの速度変動モードに対応しています．本コードでは通常のレゾルベント解析の他に，乱択アルゴリズムに基づく高速なレゾルベント解析の計算法 ([randomized resolvent analysis](https://arxiv.org/abs/1902.01458)) も実装しています．

![Figure3-2](https://user-images.githubusercontent.com/47338366/76111281-d7e30780-5f94-11ea-9f57-79398097fc97.png)


## 依存モジュール
* Python 3.*
* NumPy 
* SciPy
* [Ofpp](https://github.com/dayigu/ofpp) - OpenFOAMデータの読み込みに使用しています.
* [mpi4py](https://mpi4py.readthedocs.io/en/stable/) - オプション (推奨).

## 使用法
OpenFOAMを用いて圧縮性流れ場の計算を実施したと仮定して，本コードの使用法を説明します．なお，本コードは現在のところ，六面体格子の結果しかサポートしておりません．本コードでは流れ場の変数は次の式で無次元化されていると仮定しています．

<img src="https://latex.codecogs.com/gif.latex?x&space;=&space;\frac{\widetilde{x}}{L},&space;\:&space;y&space;=&space;\frac{\widetilde{y}}{L},&space;\:&space;z&space;=&space;\frac{\widetilde{z}}{L}">, 

<img src="https://latex.codecogs.com/gif.latex?\rho&space;=&space;\frac{\widetilde{\rho}}{\rho_\infty},&space;\:&space;u&space;=&space;\frac{\widetilde{u}}{a_\infty},&space;\:&space;T&space;=&space;\frac{\widetilde{T}}{T_\infty}">.

ここで, x, y, z は位置座標，<img src="https://latex.codecogs.com/gif.latex?\rho"> は密度，<img src="https://latex.codecogs.com/gif.latex?u"> は速度ベクトル，<img src="https://latex.codecogs.com/gif.latex?T"> は温度， <img src="https://latex.codecogs.com/gif.latex?a">は音速です．添字 <img src="https://latex.codecogs.com/gif.latex?\infty"> は一様流における値を表しています．OpenFOAMにおいて上記の無次元化を実施する場合は，等圧比熱<img src="https://latex.codecogs.com/gif.latex?c_p"> とモル質量 <img src="https://latex.codecogs.com/gif.latex?m"> を <img src="https://latex.codecogs.com/gif.latex?c_p=2.5">, <img src="https://latex.codecogs.com/gif.latex?m=11640.3"> と設定するのが最も簡単と思われます ([こちら](CylinderFlow/constant/thermophysicalProperties)の設定ファイルをご覧下さい). 

LNSオペレータを構築するには，時間平均流れ場 (`CylinderFlow/1000/*Mean`) をあらかじめ計算しておく必要があります ([例](CylinderFlow/1000)). また，本コードではセル中心座標 (`CylinderFlow/1000/C`) とセル体積 (`CylinderFlow/1000/V`) をLNSオペレータの構築に使用しています．これらのデータはOpenFOAMのケースディレクトリ内で以下のコマンドを実行することで取得できます．
```
postProcess -funcs writeCellCentres 
postProcess -funcs writeCellVolumes
```

本コードではOpenFOAMデータを読み込む際に，外部モジュールとしてOfppを使用しています．このモジュールは，(私の環境では) [controlDict](CylinderFlow/system/controlDict)内の `writeFormat` オプションを `ascii` に設定した時に最も期待通りに動作しました．Ofppモジュール内でエラーに遭遇した場合は，`writeFormat`オプションを`binary`に設定した上で`foamFormatConvert`を用いてデータをバイナリ形式に変換し，その後再びASCII形式に変換してみて下さい．


### オペレータの生成
[GenerateOperator.py](GenerateOperator.py) または [GenerateOperatorMPI.py](GenerateOperatorMPI.py) (並列計算用) を用いてオペレータを生成できます．
```
python3 GenerateOperator.py -f [Parameter file name] -p [Parameter name]
mpiexec -np [Number of thread] python3 GenerateOperatorMPI.py -f [Parameter file name] -p [Parameter name]
```

パラメータファイル ([Parameters.dat](Parameters.dat)) の設定れいを以下に示します．`CaseDir`と` TimeDir` はそれぞれシミュレーションのケースディレクトリとシミュレーションデータを含む時刻ディレクトリを設定します．
```
[GenerateOperator]
CaseDir = CylinderFlow/
TimeDir = 1000/
Operator = CylinderFlow/matL.npz
Viscosity = 1.333333e-3
PrandtlNumber = 0.7
```


### Stability or Resolvent Analysis
[CalcMode.py](CalcMode.py) または [CalcModeMPI.py](CalcModeMPI.py) を用いて安定性解析やレゾルベント解析を実行できます． 以下のコマンドを実行すると可視化用のTecPlot ASCII形式ファイル (.dat) とデータ保存用のバイナリファイル (.pickle) が出力されます．安定性解析を実施した場合は固有値，レゾルベント解析の実施時にはゲインがテキストファイルとして出力されます．また，レゾルベント解析の結果得られたゲイン分布のプロットが出力されます．ゲイン分布はその値と購買の情報を用いて3次内挿されます ([参考文献1](https://web.stanford.edu/group/ctr/Summer/SP14/08_Transition_and_turbulence/11_fosas.pdf), [参考文献2](https://spiral.imperial.ac.uk/handle/10044/1/72876))．
```
python3 CalcMode.py -f [Parameter file name] -p [Parameter name]
mpiexec -np [Number of thread] python3 CalcModeMPI.py -f [Parameter file name] -p [Parameter name]
```

このツールでは固有値問題を解く際に，shift-invertモードを使用して複素平面場の固有値を求めます．Shift-invertモードを使用する際に必要な複素数値の組みはパラメータ (以下の例では`Sigma*`) によって制御され，それぞれの複素数に対して固有値問題が解かれます．`ModeNum`は各複素数に対して計算される固有値および固有ベクトルの数を設定します．以下の例では，計算を実行すると104個 (= 4 * 26) のモードが出力されます．また，`Which`は固有値計算の優先順位を制御します (詳細は[こちら](https://docs.scipy.org/doc/scipy/reference/tutorial/arpack.html)をご覧下さい)．
```
[Stability]
Mode = Stability
CaseDir = CylinderFlow/
TimeDir = 1000/
Operator = CylinderFlow/matL.npz
SaveName = CylinderFlow/stability
ModeNum = 4
SigmaRealStart = 0.01
SigmaRealEnd = 0.01
SigmaRealNum = 1
SigmaImagStart = 0.0
SigmaImagEnd = 1.0
SigmaImagNum = 26
Which = LM
```

レゾルベン解析に際しては，計算を高速化するために [randomized method](https://arxiv.org/abs/1902.01458) の使用をおすすめします．`Omega`に関するパラメータはレゾルベント解析における周波数スイープの設定に使用されます．`Alpha`の値はレゾルベントモードの時間割引きに関係しています (詳細は[こちらの論文](https://doi.org/10.1017/jfm.2019.163)をご覧ください)．また，`ResolventMode` はどのモードを保存するかを設定しており， `Both`, `Response`, `Forcing` が使用できます.
```
[Resolvent]
Mode = RandomizedResolvent
CaseDir = CylinderFlow/
TimeDir = 1000/
Operator = CylinderFlow/matL.npz
SaveName = CylinderFlow/resolvent
ModeNum = 3
OmegaStart = 0.0
OmegaEnd = 1.0
OmegaNum = 101
AlphaStart = 0.0
AlphaEnd = 0.0
AlphaNum = 1
ResolventMode = Both
```

## 作者

* **小島 良実** - [niktFluid](https://github.com/niktFluid)


## ライセンス

This project is licensed under the MIT License - see the [LICENSE.md](LICENSE) file for details


## 謝辞

本コードはUCLAの[平研究室](http://www.seas.ucla.edu/fluidflow/index.html)訪問時に作成し，メンバーには様々な議論やアドバイスを頂きました．ここに記して感謝致します．
