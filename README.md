# Homework4-Gouraud-Shading-
Implementing graphic pipeline and rendering 3D objects by Gouraud Shading.

①開発期間： おおよそ二か月（もし2Dの部分も含めると、二か月になります）

②開発メンバー数： 一人

③使用ツール: Visual Studio, OpenGL 

③動作環境：今の運行環境はVisual Studio 15.0でコンパイルしています。また、OpenGL環境を構築する
　　　　　　必要があります。 （もしないなら、以下ウェブサイトの通りに、特定のフォルダーに入れてください。）
　　　　　　
　　　　　　リンク：http://hajimeyo-opengl.sakura.ne.jp/

④使用方法：コンパイルすると、プログラムが実行できます。もし他のモデルをチェックしたい場合、
	　　main関数load_file( "Hw4A.in" );にHw4A,　Hw4B, Hw4C,　
　　　　　　Hw4D, Hw4Eから、一つを選んで変えてください。

⑤アピールポイントと工夫した点：
OpenGL APIに頼らなく（点を描く関数以外）、3Dグラフィックスを表現するため必要なアルゴリズム、
例えば, Transform, Clipping, Hidden Surface Removal, Z-Buffer Algorithm,Gouraud Shadingなど、全部自分で最初からプログラミング
しました。
 　　　
