#include <string>
#include <memory>
#include <map>
#include <unordered_map>
#include <boost/optional.hpp>
using namespace std;

/// 値の型
typedef double val_t;
/// optional 値の型
typedef boost::optional<val_t> opt_val_t;

/// 単項演算子。
struct u_op {
	unsigned char token;
	unsigned char cost;
	function<opt_val_t(val_t)> f;
	u_op(unsigned char token, unsigned char cost, function<opt_val_t(val_t)> f) {
		this->token = token; this->cost = cost; this->f = f;
	}
};

/// 二項演算子。
struct b_op {
	unsigned char token;
	unsigned char cost;
	function<opt_val_t(val_t, val_t)> f;
	b_op(unsigned char token, unsigned char cost, function<opt_val_t(val_t, val_t)> f) {
		this->token = token; this->cost = cost; this->f = f;
	}
};

/// 演算子の特性
struct op_c {
	/// 単項演算ならtrue、二項演算ならfalse
	bool unary;

	string prefix;
	string infix;
	string suffix;

	/// 優先度
	int precedence;
	/// 左結合性（b_opのみ）
	bool left_associative;
	/// 右結合性（b_opのみ）
	bool right_associative;

	/// 括弧が必要（u_opのみ）
	bool needs_paren;
	/// オペランド１に括弧が必要（b_opのみ）
	bool needs_paren_op1;
	/// オペランド２に括弧が必要（b_opのみ）
	bool needs_paren_op2;

	/// 単項演算子のためのコンストラクタ
	op_c(string prefix, string suffix, int precedence, bool needs_paren) {
		this->unary = true;
		this->prefix = prefix;
		this->suffix = suffix;
		this->precedence = precedence;
		this->needs_paren = needs_paren;
	}

	/// 二項演算子のためのコンストラクタ
	op_c(string prefix, string infix, string suffix,
		int precedence, bool left_associative, bool right_associative,
		bool needs_paren_op1, bool needs_paren_op2) {
			this->unary = false;
			this->prefix = prefix; this->infix = infix; this->suffix = suffix;
			this->precedence = precedence; this->left_associative = left_associative; this->right_associative = right_associative;
			this->needs_paren_op1 = needs_paren_op1; this->needs_paren_op2 = needs_paren_op2;
	}

	op_c() {}
};

/// 式表現
typedef unsigned char exp_rep_t[];
/// 式表現のunique_ptr
typedef unique_ptr<exp_rep_t> exp_rep_p;

/*
式表現の形式：
第０バイト：トークンの個数
あとは4bitずつのトークン列
奇数個の場合は上位4bitが最後のトークン、残りは0。
*/

/// 式表現の指定位置のトークンを取得
unsigned char token_at(exp_rep_t e, int i);
/// 式表現の長さを取得
size_t explen(exp_rep_t e);
/// 式表現のトークン列
string expstr(exp_rep_t e);
/// 式表現のTeX文字列
string exptex(exp_rep_t e, string seed);

/// 数値の表現だけを作る
exp_rep_p n_combine(size_t length);
/// 単項演算で式表現を延長
exp_rep_p u_combine(exp_rep_t e, const unsigned char t);
/// 二項演算で式表現を結合
exp_rep_p b_combine(exp_rep_t e1, exp_rep_t e2, const unsigned char t);

/// 式。
struct exp_t {
	/// 表現
	exp_rep_p rep;
	/// コスト
	unsigned char cost;
};


/// 等式の集合
typedef unordered_map<val_t, exp_t> equations;

/// 与えられた文字列から生成される等式の集合を計算する
equations&& think(string numbers);
/// 単項演算で拡大 seedから（poolに含まれていないものを）opで拡大
equations u_extend(equations&& seed, equations&& pool, u_op op);
/// 二項演算で拡大
equations b_extend(equations&& lhs, equations&& rhs, b_op op);


/// 整数っぽい値なら整数化した値を、そうでないならばnoneを
opt_val_t maybe_int(val_t);
/// 整数っぽい値なら整数化した値を、そうでないならば元の値を、noneならnoneを
opt_val_t round_if_needed(opt_val_t);
/// 値のカットオフ（大きすぎる途中結果の排除）
bool cutoff(val_t);
