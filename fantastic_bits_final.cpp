#pragma GCC optimize("O3") 
#pragma GCC optimize("inline")
#pragma GCC optimize("omit-frame-pointer")
#include <ctime>

#define DEPTH 3
#define INIT_RANDOM_NUM 400
#define PI_CONST 3.141592
#define POOL_FIT_FORMULA current_fit

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <utility> //pair
#include <map>
#include <cmath> //sqrt

typedef std::pair<int, int> vec2;
typedef std::pair<double, double> vec2f;

struct entity{
	int uid, pid, state;
	enum entity_type{ wizard, snaffle, bludger, immobile } type;
	double x, y, vx, vy;
	int grabbed_snaffle_uid;
};

struct bounds_collision_data{ double time; int side; };

struct move{
	enum move_type{
		MOVE, THROW, OBLIVIATE, PETRIFICUS, ACCIO, FLIPENDO
	} type;
	int x, y, thrust, t_uid;
};

typedef std::pair<move, move> moves;

class gameState{
	struct spell{
		enum spell_type{
			OBLIVIATE, PETRIFICUS, ACCIO, FLIPENDO
		} type;
		int user, target;
		int timer;
	};

	std::vector<spell> obliviates_vec;
	std::vector<spell> petrificus_vec;
	std::vector<spell> accio_vec;
	std::vector<spell> flipendo_vec;
public:
	gameState();
	void init();
	void play(const moves& move_player_A, const moves& move_player_B);

	int score_A, score_B;
	int mana_A, mana_B;
	std::vector<entity> entities;
	std::map<int, entity*> uid_map;

	vec2f closest(const vec2f& p, const vec2f& l1, const vec2f& l2);

	double collision(const entity& a, const vec2f& av, const entity& b, const vec2f& bv);
	void resolve_collision(entity*, entity*);

	bounds_collision_data bounds_collision(entity);
	void resolve_bounds_collision(bounds_collision_data, entity*);
};



double get_radius(const entity&);
double get_mass(const entity&);

double euclid_dist2(const vec2f& V);

vec2f get_direction(const float& angle, const vec2f& pos);

class dummy{
public:
	void game_start(gameState*, int);
	void game_end();
	moves take_turn();
};

class simple_AI{
	int my_goal_x, my_goal_y;
	gameState *state;
	int id;
public:
	void game_start(gameState*, int);
	void game_end();
	moves take_turn();
};

struct gene {
	float throw_arc; //rad, angles for consecutive throws
	float move_arcs[DEPTH]; //rad, angles to move towards
	int accio_target, flip_target; //uid mod number of potential targets
	int accio_turn, flip_turn;
};

struct chromosome {
	gene A, B;
};

class Genetic_A{
	chromosome genes_random();
	float fitness(const chromosome&);

	moves extract_moves(const chromosome& T, int turn, gameState&); //assumes entity is at 0,0

	vec2f my_goal;
	gameState *state;
	int my_id;
	int debug;
public:
	void game_start(gameState*, int);
	void game_end();
	moves take_turn();
};

std::pair<chromosome, chromosome> roulette_pair(const std::vector<std::pair<float, chromosome>>& pool);
chromosome gene_crossover(const chromosome& A, const chromosome& B);
chromosome gene_mutate(const chromosome& A);


double euclid_dist2(const vec2f& V){
	return (V.first*V.first + V.second*V.second);
}

vec2f resize(const vec2f& V, const double& scalar){
	auto d = sqrt(euclid_dist2(V));
	return vec2f((scalar*V.first)/d, (scalar*V.second)/d);
}

gameState::gameState(){
	score_A = 0;
	score_B = 0;
}

void gameState::init(){
	for(auto &iv : this->entities){
		this->uid_map[iv.uid] = &iv;
	}
}

void gameState::play(const moves& move_player_A, const moves& move_player_B) {

	for(auto &iv : this->entities){
		// apply wizard move
		if(iv.type == entity::entity_type::wizard){
			move wiz_move;
			if(iv.pid == 0){
				if(iv.uid == 0 || iv.uid == 2) wiz_move = move_player_A.first; // first wizard A
				else wiz_move = move_player_A.second; // second wizard A
			}
			else{
				if(iv.uid == 0 || iv.uid == 2) wiz_move = move_player_B.first; // first wizard B
				else wiz_move = move_player_B.second; // second wizard B
			}

			if(wiz_move.type == move::move_type::MOVE && wiz_move.thrust > 0){
				vec2f wiz_move_vec(wiz_move.x - iv.x, wiz_move.y - iv.y);
				wiz_move_vec = resize(wiz_move_vec, wiz_move.thrust * 1.0); // *1.0 represents the mass of the wizard

				iv.vx = uid_map[iv.uid]->vx + wiz_move_vec.first;
				iv.vy = uid_map[iv.uid]->vy + wiz_move_vec.second;
			}
			
			else if(wiz_move.type == move::move_type::THROW){
				//if(iv.state == 0)
					//incorrect input, end the game and such
				vec2f throw_vec(wiz_move.x - iv.x, wiz_move.y - iv.y);
				throw_vec = resize(throw_vec, wiz_move.thrust * 2.0); // division by 0.5 of the snaffle mass
				uid_map[iv.grabbed_snaffle_uid]->vx += throw_vec.first;
				uid_map[iv.grabbed_snaffle_uid]->vy += throw_vec.second;
				uid_map[iv.grabbed_snaffle_uid]->state = 0;
			}
			
			else if(wiz_move.type == move::move_type::OBLIVIATE){
				if(iv.pid == 0){
					if(mana_A < 5) continue;
					mana_A -= 5;
				}
				else if(iv.pid == 1){
					if(mana_B < 5) continue;
					mana_B -= 5;
				}
				this->obliviates_vec.push_back(spell{spell::OBLIVIATE, iv.pid, wiz_move.t_uid, 4});
			}

			else if(wiz_move.type == move::move_type::PETRIFICUS){
				if(iv.pid == 0){
					if(mana_A < 10) continue;
					mana_A -= 10;
				}
				else if(iv.pid == 1){
					if(mana_B < 10) continue;
					mana_B -= 10;
				}
				this->petrificus_vec.push_back(spell{spell::PETRIFICUS, iv.uid, wiz_move.t_uid, 3});
			}

			else if(wiz_move.type == move::move_type::ACCIO){
				if(iv.pid == 0){
					if(mana_A < 20) continue;
					mana_A -= 20;
				}
				else if(iv.pid == 1){
					if(mana_B < 20) continue;
					mana_B -= 20;
				}
				this->accio_vec.push_back(spell{spell::ACCIO, iv.uid, wiz_move.t_uid, 7});
			}

			else if(wiz_move.type == move::move_type::FLIPENDO){
				if(iv.pid == 0){
					if(mana_A < 20) continue;
					mana_A -= 20;
				}
				else if(iv.pid == 1){
					if(mana_B < 20) continue;
					mana_B -= 20;
				}
				this->flipendo_vec.push_back(spell{spell::FLIPENDO, iv.uid, wiz_move.t_uid, 4});
			}
		} 

		// bludger logic and move
		// for bludgers, state stores uid of the last wizard it collided with
		if(iv.type == entity::entity_type::bludger){
			double closest = 99999;
			entity target;
			target.uid = -1;
			for(auto kv : this->entities){ // find a target
				if(kv.type == entity::entity_type::wizard && kv.uid != iv.state){
					bool obliviated = false;
					for(auto obl_v : obliviates_vec){
						if(obl_v.target == iv.uid && kv.pid == obl_v.user && obl_v.timer < 4){
							obliviated = true;
							break;
						}
					}
					if(obliviated) continue;
					double d = sqrt(euclid_dist2(vec2f(iv.x - kv.x, iv.y - kv.y)));
					if(d < closest){
						closest = d;
						target = kv;
					}
				}
			}

			if(target.uid != -1){
				// found target, calculate change vector and apply thrust
				vec2f bludger_move_vec(target.x - iv.x, target.y - iv.y);
				bludger_move_vec = resize(bludger_move_vec, 1000.0 / 8.0); // bludger mass = 8, bludger thrust = 1000
				iv.vx = uid_map[iv.uid]->vx + bludger_move_vec.first;
				iv.vy = uid_map[iv.uid]->vy + bludger_move_vec.second;
			}
			// otherwise no valid targets, don't do anything
		}
	}

	// apply spell power
	auto ptrf_sv = petrificus_vec.begin();
	while(ptrf_sv != petrificus_vec.end()){
		ptrf_sv->timer--;
		if(ptrf_sv->timer == 1){
			uid_map[ptrf_sv->target]->vx = 0;
			uid_map[ptrf_sv->target]->vy = 0;
		}
		if(ptrf_sv->timer == 0){
			ptrf_sv = petrificus_vec.erase(ptrf_sv);
		}
		else if(ptrf_sv != petrificus_vec.end()) ptrf_sv++;
	}

	auto acc_sv = accio_vec.begin();
	while(acc_sv != accio_vec.end()){
		acc_sv->timer--;
		if(acc_sv->timer < 6){
			vec2f dir(uid_map[acc_sv->target]->x - uid_map[acc_sv->user]->x, uid_map[acc_sv->target]->y - uid_map[acc_sv->user]->y);

			double d = euclid_dist2(dir);
			if(d > 0){
				double spell_power = 3000.0 / (d / (1000.0*1000.0));
				spell_power = spell_power < 1000.0 ? spell_power : 1000.0;
				spell_power /= get_mass(*uid_map[acc_sv->target]);
				dir = resize(dir, spell_power);
				uid_map[acc_sv->target]->vx -= dir.first;
				uid_map[acc_sv->target]->vy -= dir.second;
			}
			else{
				acc_sv = accio_vec.erase(acc_sv);
				continue;
			}
		}
		if(acc_sv->timer == 0){
			acc_sv = accio_vec.erase(acc_sv);
		}
		else if(acc_sv != accio_vec.end()) acc_sv++;
	}

	auto sv = flipendo_vec.begin();
	while(sv != flipendo_vec.end()){
		sv->timer--;
		if(sv->timer < 3){
			vec2f dir(uid_map[sv->target]->x - uid_map[sv->user]->x, uid_map[sv->target]->y - uid_map[sv->user]->y);
			
			double d = euclid_dist2(dir);
			if(d > 0){
				double spell_power = 6000.0 / (d / (1000.0*1000.0));
				spell_power = spell_power < 1000.0 ? spell_power : 1000.0;
				spell_power /= get_mass(*uid_map[sv->target]);
				dir = resize(dir, spell_power);

				uid_map[sv->target]->vx += dir.first;
				uid_map[sv->target]->vy += dir.second;
			}
		}
		if(sv->timer == 0){
			sv = flipendo_vec.erase(sv);
		}
		else if(sv != flipendo_vec.end()) sv++;
	}

	auto o_sv = obliviates_vec.begin();
	while(o_sv != obliviates_vec.end()){
		o_sv->timer--;
		if(o_sv->timer == 0) o_sv = obliviates_vec.erase(o_sv);
		else if(o_sv != obliviates_vec.end()) o_sv++;
	}

	// resolve collisons
	double turn_time_left = 1.0;
	// 1. go through every pair of entities
	// 2. find the collision that will happen first
	// 3. resolve that collision, update all pods' position to that time frame
	// 4. go back to 1, repeat as long as there are collisions left
	while(1){
		bounds_collision_data tmp_col, earliest_b_col;
		double earliest_collision(2);
		entity *col_entity_A(nullptr), *col_entity_B(nullptr);
		
		for (int i = 0; i < entities.size(); ++i) {
			for (int k = i + 1; k < entities.size(); ++k) {
				double t = collision(entities[i], vec2f(entities[i].vx, entities[i].vy),
					entities[k], vec2f(entities[k].vx, entities[k].vy));
				if (t > -0.1 && t < earliest_collision) {
					earliest_collision = t;
					col_entity_A = &entities[i];
					col_entity_B = &entities[k];
				}
			}
			tmp_col = bounds_collision(entities[i]);
			if (tmp_col.time > -0.1 && tmp_col.time < earliest_collision) {
				earliest_collision = tmp_col.time;
				col_entity_A = &entities[i];
				col_entity_B = nullptr;
				earliest_b_col = tmp_col;
			}
		}
		if(earliest_collision > turn_time_left) break;

		for(auto &gv : this->entities){
			gv.x += gv.vx*earliest_collision;
			gv.y += gv.vy*earliest_collision;
		}

		turn_time_left -= earliest_collision;

		if(col_entity_B == nullptr){
			resolve_bounds_collision(earliest_b_col, col_entity_A);
		}
		else resolve_collision(col_entity_A, col_entity_B);
	}
	// update positions
	for(auto &iv : this->entities){
		iv.x += iv.vx*turn_time_left;
		iv.y += iv.vy*turn_time_left;
	}

	for(auto &iv : this->entities){
		if(iv.type == entity::wizard && iv.state == 4){
			uid_map[iv.grabbed_snaffle_uid]->vx = iv.vx;
			uid_map[iv.grabbed_snaffle_uid]->vy = iv.vy;
			uid_map[iv.grabbed_snaffle_uid]->x = iv.x;
			uid_map[iv.grabbed_snaffle_uid]->y = iv.y;
		}
	}

	// apply friction
	for(auto &iv : this->entities){
		if(iv.type == entity::entity_type::wizard || iv.type == entity::entity_type::snaffle){
			iv.vx = iv.vx*0.75;
			iv.vy = iv.vy*0.75;
			if(iv.type == entity::entity_type::snaffle || iv.type == entity::wizard){
				if(iv.state > 0){
					iv.state--;
				}
			}
		}
		else if(iv.type == entity::entity_type::bludger){
			iv.vx = iv.vx*0.9;
			iv.vy = iv.vy*0.9;
		}

		iv.x = std::round(iv.x);
		iv.y = std::round(iv.y);
		iv.vx = std::round(iv.vx);
		iv.vy = std::round(iv.vy);
	}
	mana_A++;
	mana_B++;
}

double get_radius(const entity& a){
	if(a.type == entity::snaffle) return 150.0;
	if(a.type == entity::bludger) return 200.0;
	if(a.type == entity::immobile) return 300.0;
	return 400.0;
}

double get_mass(const entity& a){
	if(a.type == entity::snaffle) return 0.5;
	if(a.type == entity::bludger) return 8;
	if(a.type == entity::immobile) return 9999999;
	return 1;
}

bounds_collision_data gameState::bounds_collision(entity a){
	int closest_side(-1);
	double earliest_collision(999999);
	
	double x1(a.x), y1(a.y), x2(a.x + a.vx), y2(a.y + a.vy);
	
	for(int side(0); side < 4; ++side){
		double x3, y3, x4, y4;
		if(side == 0){
			x3 = y3 = y4 = 0; x4 = 16000;
		}
		else if(side == 1){
			x3 = x4 = 16000; y3 = 0; y4 = 7500;
		}
		else if(side == 2){
			x3 = 0; x4 = 16000; y3 = y4 = 7500;
		}
		else{
			x3 = x4 = y3 = 0; y4 = 7500;
		}

		// find (Px, Py), the point of intersection of the line and travel path of the circle
		double det = ((x1 - x2)*(y3 - y4) - (y1 - y2)*(x3 - x4));
		
		if(det == 0) continue;

		double Px = ((x1*y2 - y1*x2)*(x3 - x4) - (x1 - x2)*(x3*y4 - y3*x4)) / det;
		double Py = ((x1*y2 - y1*x2)*(y3 - y4) - (y1 - y2)*(x3*y4 - y3*x4)) / det;

		//check if we're going towards the intersection point
		if( !
			((side == 0 && a.vy < 0)
			|| (side == 1 && a.vx > 0)
			|| (side == 2 && a.vy > 0)
			|| (side == 3 && a.vx < 0)
			))
			continue;

		// find the angle we're hitting the wall at
		double e;
		if(side == 0 || side == 2) e = atan(a.vy / a.vx);
		// rotate by 90 degrees if checking with the vertical walls
		else if(side == 1 || side == 3) e = atan(a.vx / a.vy);
		double rad = get_radius(a);

		if(a.type == entity::snaffle && (side == 1 || side == 3) && Py > 1750 && Py < 5750){
			rad = 0;
		}

		// finally the distance from Px, Py to collision point
		double c = rad / sin(e);
		vec2f to_col(Px - x1, Py - y1); // vector from a to intersection
		to_col = resize(to_col, sqrt(euclid_dist2(to_col)) - c);
		vec2f to_col_sec(Px - x1, Py - y1);
		to_col_sec = resize(to_col_sec, sqrt(euclid_dist2(to_col)) + c + c);
		double distance_to_col = sqrt(euclid_dist2(to_col));
		double distance_to_col_sec = sqrt(euclid_dist2(to_col_sec));
		double travel_distance = sqrt(euclid_dist2(vec2f(a.vx, a.vy)));

		double time = distance_to_col < distance_to_col_sec ? distance_to_col : distance_to_col_sec;
		time /= travel_distance;

		if(time < earliest_collision){
			earliest_collision = time;
			closest_side = side;
		}
	}
	return{earliest_collision, closest_side};
};

void gameState::resolve_bounds_collision(bounds_collision_data col, entity *a){
	if(col.side == 0 || col.side == 2){
		a->vy = -a->vy;
	}
	else{
		a->vx = -a->vx;

		//scoring
		if(a->type == entity::snaffle){
			if(a->y > 1750 && a->y < 5750){
				for(int i(0); i < entities.size(); ++i){
					if(entities[i].uid == a->uid){
						entities.erase(entities.begin() + i);
						break;
					}
				}
				//assign point
				if(col.side == 1) this->score_A++;
				else this->score_B++;
			}
		}

	}
}

/*
* the following collision prediction and resolution code was written by Magus for the Coders Strike Back contest
* files.magusgeek.com/csb/csb_en.html
*/

vec2f gameState::closest(const vec2f& p, const vec2f& l1, const vec2f& l2) {
	double da = l2.second - l1.second;
	double db = l1.first - l2.first;
	double c1 = da*l1.first + db*l1.second;
	double c2 = -db*p.first + da*p.second;
	double det = da*da + db*db;
	float cx = 0, cy = 0;
	if(det != 0){
		cx = (da*c1 - db*c2) / det;
		cy = (da*c2 + db*c1) / det;
	}
	else{
		cx = p.first;
		cy = p.second;
	}
	return vec2f(cx, cy);
}

double gameState::collision(const entity& a, const vec2f& av, const entity& b, const vec2f& bv){
	// ignore snaffles if they were grabbed this turn
	if(a.type == entity::snaffle && a.state == 1){
		return -1;
	}
	else if(b.type == entity::snaffle && b.state == 1){
		return -1;
	}
	// ignore snaffle-wizard collisions if wizard grabbed a snaffle recently
	else if(a.type == entity::wizard && b.type == entity::snaffle && a.state > 1){
		return -1;
	}
	else if(b.type == entity::wizard && a.type == entity::snaffle && b.state > 1){
		return -1;
	}

	double dist2 = euclid_dist2(vec2f(a.x - b.x, a.y - b.y));
	double rads2 = get_radius(a) + get_radius(b);

	if(a.type == entity::wizard && b.type == entity::snaffle) rads2 = 399.0;
	else if(b.type == entity::wizard && a.type == entity::snaffle) rads2 = 399.0;

	rads2 *= rads2;

	vec2f a_(a.x - b.x, a.y - b.y);
	vec2f av_(av.first - bv.first, av.second - bv.second);
	vec2f b_(0, 0);
	
	vec2f p = closest(b_, a_, vec2f(av_.first + a_.first, av_.second + a_.second));
	double pdist2 = euclid_dist2(p);
	double a_pdist2 = euclid_dist2(vec2f(a_.first - p.first, a_.second - p.second));

	if(pdist2 < rads2){
		double length = sqrt(euclid_dist2(av_));
		double backdist = sqrt(rads2 - pdist2);
		p.first = p.first - ((backdist*av_.first) / length);
		p.second = p.second - ((backdist*av_.second) / length);
		
		double pdist2 = euclid_dist2(vec2f(a_.first - p.first, a_.second - p.second));
		if(pdist2 > a_pdist2) return -1;
		double pdist = sqrt(pdist2);
		if(pdist > length) return -1;

		return pdist / length; //time to collision
	}
	return -1;
}

void gameState::resolve_collision(entity *A, entity *B){
	if(A->type == entity::snaffle && B->type == entity::wizard){
		A->vx = B->vx, A->vy = B->vy, A->x = B->x, A->y = B->y;
		B->state = 4;
		A->state = 1;
		B->grabbed_snaffle_uid = A->uid;
		return;
	}
	else if(A->type == entity::wizard && B->type == entity::snaffle){
		B->vx = A->vx, B->vy = A->vy, B->x = A->x, B->y = A->y;
		A->state = 4;
		B->state = 1;
		A->grabbed_snaffle_uid = B->uid;
		return;
	}
	// if a bludger hit someone, save their UID in the bludger's state
	if(A->type == entity::bludger && B->type == entity::wizard){
		A->state = B->uid;
	}
	else if(B->type == entity::bludger && A->type == entity::wizard){
		B->state = A->uid;
	}
	double mA = get_mass(*A);
	double mB = get_mass(*B);

	double mcoeff = (mA + mB) / (mA * mB);

	vec2f n(A->x - B->x, A->y - B->y);
	double nsq(euclid_dist2(n));

	double dvx = A->vx - B->vx;
	double dvy = A->vy - B->vy;

	double product = n.first*dvx + n.second*dvy;
	double fx = (n.first*product) / (nsq * mcoeff);
	double fy = (n.second*product) / (nsq * mcoeff);

	double impulse = sqrt(fx*fx + fy*fy);

	A->vx -= fx / mA;
	A->vy -= fy / mA;
	B->vx += fx / mB;
	B->vy += fy / mB;
	
	if(impulse < 100.0){
		fx *= 100.0 / impulse;
		fy *= 100.0 / impulse;
	}
	A->vx -= fx / mA;
	A->vy -= fy / mA;
	B->vx += fx / mB;
	B->vy += fy / mB;

	if(A->type == entity::immobile){
		A->vx = 0;
		A->vy = 0;
	}
	else if(B->type == entity::immobile){
		B->vx = 0;
		B->vy = 0;
	}
}

vec2f get_direction(const float& angle, const vec2f& pos) {
	return vec2f(1000 * cos(angle) + pos.first, 1000 * sin(angle) + pos.second);
}

moves simple_AI::take_turn(){
	bool first_move = 0;
	move a, b;

	int mana;
	if(this->id == 0) mana = state->mana_A;
	else mana = state->mana_B;

	int busy_target(-1), potential_accio_goalkeep(-1);
	int farthest_snaffle(0);
	for (auto kv : state->entities) {
		if (kv.type == entity::snaffle) {
			int d = abs(my_goal_x - kv.x);
			if (d > farthest_snaffle && d > 10000) {
				farthest_snaffle = d;
				potential_accio_goalkeep = kv.uid;
			}
		}
	}

	for(auto iv : state->entities){
		if(iv.type == entity::wizard && iv.pid == this->id){
		move curr;

		if(iv.state == 3){
			curr.type = move::THROW;
			curr.thrust = 500;
			curr.x = this->my_goal_x;
			curr.y = this->my_goal_y;
		} //throwing

		else if(mana > 19 && potential_accio_goalkeep != -1){
			curr.type = move::ACCIO;
			mana -= 20;
			curr.t_uid = potential_accio_goalkeep;
		}
		else{
			int closest_snaffle = 0, d = 99999999;
			for(auto kv : state->entities){
				if(kv.type == entity::snaffle){
					if (kv.uid == busy_target) continue;
					int d_ = (int)euclid_dist2(vec2f(iv.x - kv.x, iv.y - kv.y));
					if(d_ < d){
						closest_snaffle = kv.uid;
						d = d_;
					}
				}
			} // find closest snaffle loop
			
			curr.type = move::MOVE;
			curr.x = (int)state->uid_map[closest_snaffle]->x;
			curr.y = (int)state->uid_map[closest_snaffle]->y;
			curr.thrust = 150;
			if (busy_target == -1) busy_target = closest_snaffle;
		} // moving


		if(!first_move){
			a = curr, first_move = 1;
		}
		else b = curr;
	}
	}
	return moves(a, b);
}

void simple_AI::game_start(gameState* state_reference, int my_id){
	state = state_reference;
	id = my_id;

	if(id == 0){
		my_goal_x = 16000;
		my_goal_y = 3750;
	}
	else{
		my_goal_x = 0;
		my_goal_y = 3750;
	}
}

void simple_AI::game_end(){}

void dummy::game_end(){}
void dummy::game_start(gameState*, int){}
moves dummy::take_turn(){
	move a{move::MOVE, 8000, 3750, 0, 0};
	return moves(a, a);
}


chromosome gene_crossover(const chromosome& A, const chromosome& B) {
	if(rand() % 2){
		chromosome T = A;
		T.B = B.B;
		return T;
	}
	else {
		chromosome T = B;
		T.B = A.B;
		return T;
	}
}


chromosome gene_mutate(const chromosome& A) {
	return chromosome(A);
}

chromosome Genetic_A::genes_random()
{
	gene A, B;
	for (int i(0); i < DEPTH; i++) {
		A.move_arcs[0] = static_cast <float> (rand()) / static_cast <float> (RAND_MAX / (2 * PI_CONST)) - PI_CONST;
		B.move_arcs[0] = static_cast <float> (rand()) / static_cast <float> (RAND_MAX / (2 * PI_CONST)) - PI_CONST;
	}
	A.throw_arc = static_cast <float> (rand()) / static_cast <float> (RAND_MAX / (2 * PI_CONST)) - PI_CONST;
	A.accio_target = rand() % 10000;
	A.accio_turn = rand() % DEPTH;
	A.flip_target = rand() % 10000;
	A.flip_turn = rand() % DEPTH;
	B.throw_arc = static_cast <float> (rand()) / static_cast <float> (RAND_MAX / (2 * PI_CONST)) - PI_CONST;
	B.accio_target = rand() % 10000;
	B.accio_turn = rand() % DEPTH;
	B.flip_target = rand() % 10000;
	B.flip_turn = rand() % DEPTH;
	return{ A, B };
}

std::pair<chromosome, chromosome> roulette_pair(const std::vector<std::pair<float, chromosome>>& pool) {
	float total(0);

	for (auto i(pool.begin()); i != pool.end(); ++i) {
		total += i->first;
	}

	float r1 = total * ((float)rand() / RAND_MAX);
	float r2 = total * ((float)rand() / RAND_MAX);
	float k = 0.;

	chromosome g1, g2;
	for (auto i(pool.begin()); i != pool.end(); ++i) {
		if (k + i->first >= r1) g1 = i->second;
		if (k + i->first >= r2) g2 = i->second;
		if (k + i->first >= r2 && k + i->first >= r1) {
			return std::make_pair(g1, g2);
		}

		k += i->first;

	}
	throw "ummmmmmm";
}


void Genetic_A::game_end() {}

void Genetic_A::game_start(gameState* state_ref, int id) {
	state = state_ref, my_id = id;
	if (id == 0) {
		my_goal = vec2f(16000, 3750);
	}
	else my_goal = vec2f(0, 3750);
}

moves Genetic_A::take_turn() {
	clock_t start = clock();
	this->debug = 0;

	std::vector<std::pair<float, chromosome>> population;
	chromosome current_g, best_g; float current_fit, best_fit(-9999);

	for (int i(0); i < INIT_RANDOM_NUM; ++i) {
		current_g = genes_random();
		current_fit = fitness(current_g);
		if (current_fit > best_fit) {
			best_fit = current_fit;
			best_g = current_g;
		}
		population.push_back(std::make_pair(current_fit, current_g));
	}

	while (double(clock() - start) / double(CLOCKS_PER_SEC) < 0.09) {

		std::pair<chromosome, chromosome> Gcp = roulette_pair(population);
		current_g = gene_mutate(gene_crossover(Gcp.first, Gcp.second));
		current_fit = fitness(current_g);

		if (current_fit > best_fit) {
			best_fit = current_fit, best_g = current_g;
		}

		population.push_back(std::make_pair(POOL_FIT_FORMULA, current_g));
	}

	this->debug = 0;

	fitness(best_g);

	if (debug) {
		auto clone = best_g;
		clone.B.throw_arc = -0.6;
		fitness(clone);
		std::cerr << population.size() << " " << best_fit << std::endl;
		std::cerr << "athrow: " << extract_moves(clone, 0, *state).second.x << "x y"
			<< extract_moves(clone, 0, *state).second.y << std::endl;
	}
	std::cerr << best_fit << ", best of " << population.size() << std::endl;

	return extract_moves(best_g, 0, *state);
}

moves Genetic_A::extract_moves(const chromosome& T, int turn, gameState& src) {
	move mov_A, mov_B;
	int mana;
	if (this->my_id == 0) mana = state->mana_A;
	else mana = state->mana_B;

	entity *A_dat, *B_dat;
	if (my_id == 0) {
		A_dat = src.uid_map[0];
		B_dat = src.uid_map[1];
	}
	else {
		A_dat = src.uid_map[2];
		B_dat = src.uid_map[3];
	}

	int potential_targets(0);
	for (auto iv : src.entities) {
		if (iv.type == entity::snaffle) {
			potential_targets++;
		}
	}

	if (A_dat->state == 3) {
		mov_A.type = move::THROW;
		mov_A.thrust = 500;
		vec2f dir = get_direction(T.A.throw_arc, vec2f(A_dat->x, A_dat->y));
		mov_A.x = dir.first;
		mov_A.y = dir.second;
	}
	else if (T.A.flip_turn == turn && mana >= 20 && potential_targets > 0) {
		mov_A.type = move::FLIPENDO;
		int target(0);
		for (auto iv : src.entities) {
			if (iv.type == entity::snaffle) {
				target++;
			}
		}
		target = T.A.flip_target % target;
		for (auto iv : src.entities) {
			if (iv.type == entity::snaffle) {
				if (!target) {
					target = iv.uid;
					break;
				}
				target--;
			}
		}
		mov_A.t_uid = target;
	}
	else if (T.A.accio_turn == turn && mana >= 20 && potential_targets > 0) {
		mov_A.type = move::ACCIO;
		int target(0);
		for (auto iv : src.entities) {
			if (iv.type == entity::snaffle) {
				target++;
			}
		}
		target = T.A.accio_target % target;
		for (auto iv : src.entities) {
			if (iv.type == entity::snaffle) {
				if (!target) {
					target = iv.uid;
					break;
				}
				target--;
			}
		}
		mov_A.t_uid = target;
	}
	else {
		mov_A.type = move::MOVE;
		mov_A.thrust = 150;
		vec2f dir = get_direction(T.A.move_arcs[turn], vec2f(A_dat->x, A_dat->y));
		mov_A.x = dir.first;
		mov_A.y = dir.second;
	}

	if (B_dat->state == 3) {
		mov_B.type = move::THROW;
		mov_B.thrust = 500;
		vec2f dir = get_direction(T.B.throw_arc, vec2f(B_dat->x, B_dat->y));
		mov_B.x = dir.first;
		mov_B.y = dir.second;
	}
	else if (T.B.flip_turn == turn && mana >= 20 && potential_targets > 0) {
		mov_B.type = move::FLIPENDO;
		int target(0);
		for (auto iv : src.entities) {
			if (iv.type == entity::snaffle) {
				target++;
			}
		}
		target = T.B.flip_target % target;
		for (auto iv : src.entities) {
			if (iv.type == entity::snaffle) {
				if (!target) {
					target = iv.uid;
					break;
				}
				target--;
			}
		}
		mov_B.t_uid = target;
	}
	else if (T.B.accio_turn == turn && mana >= 20 && potential_targets > 0) {
		mov_B.type = move::ACCIO;
		int target(0);
		for (auto iv : src.entities) {
			if (iv.type == entity::snaffle) {
				target++;
			}
		}
		target = T.B.accio_target % target;
		for (auto iv : src.entities) {
			if (iv.type == entity::snaffle) {
				if (!target) {
					target = iv.uid;
					break;
				}
				target--;
			}
		}
		mov_B.t_uid = target;
	}
	else {
		mov_B.type = move::MOVE;
		mov_B.thrust = 150;
		vec2f dir = get_direction(T.B.move_arcs[turn], vec2f(B_dat->x, B_dat->y));
		mov_B.x = dir.first;
		mov_B.y = dir.second;
	}

	return moves(mov_A, mov_B);
}

float Genetic_A::fitness(const chromosome& specimen) {
	gameState fake = *(this->state);
	fake.init();
	//fake.mana_A = 0;
	//fake.mana_B = 0;
	
	simple_AI dumb;

	dumb.game_start(&fake, my_id == 0 ? 1 : 0);

	int move_i(0);

	for (int i(0); i < DEPTH; ++i) {
		moves my_moves = extract_moves(specimen, i, fake);

		if (my_id == 0) fake.play(my_moves, dumb.take_turn());
		else fake.play(dumb.take_turn(), my_moves);
	}

	//fitness = a*score + b*fn(distances to snaffles) + c*gn(distances of snaffles to goal)
	float a(1000), b(1000), c(1);
	float fitness_value = 0;

	if (my_id == 0) fitness_value += a*(fake.score_A);
	else fitness_value += a*(fake.score_B);

	entity *A_dat, *B_dat;
	if (my_id == 0) {
		A_dat = fake.uid_map[0];
		B_dat = fake.uid_map[1];
	}
	else {
		A_dat = fake.uid_map[2];
		B_dat = fake.uid_map[3];
	}
    
    int A_snaffle(-2), B_snaffle(-2);
	double A_snaffle_val(0), B_snaffle_val(0);
	double A_best_val(0), B_best_val(0);

	for (auto iv : fake.entities) {
		if (iv.type == entity::snaffle) {
			float d_goal = sqrt(euclid_dist2(vec2f(iv.x - my_goal.first, iv.y - my_goal.second)));
			fitness_value += b * (40 * 1000) / d_goal;

			float snaffle_value = 0;
			for (auto kv : fake.entities) {
				if (kv.type == entity::snaffle && kv.uid != iv.uid) {
					//clusters
					float nearby_snaffle_d = sqrt(euclid_dist2(vec2f(iv.x - kv.x, iv.y - kv.y)));
					snaffle_value += sqrt(nearby_snaffle_d);
				}
			}
			snaffle_value = 10000.0 / snaffle_value;
			snaffle_value += d_goal / 100.0;

			A_snaffle_val = sqrt(euclid_dist2(vec2f(iv.x - A_dat->x, iv.y - A_dat->y)));
			if (A_snaffle_val < 400) A_snaffle_val = 400;
			A_snaffle_val = snaffle_value + (c * 40 * 1000 / A_snaffle_val);
			if (A_best_val < A_snaffle_val) {
				A_best_val = A_snaffle_val;
				A_snaffle = iv.uid;
			}
			B_snaffle_val = sqrt(euclid_dist2(vec2f(iv.x - B_dat->x, iv.y - B_dat->y)));
			if (B_snaffle_val < 400) B_snaffle_val = 400;
			B_snaffle_val = snaffle_value + (c * 40 * 1000 / B_snaffle_val);
			if (B_best_val < B_snaffle_val) {
				B_best_val = B_snaffle_val;
				B_snaffle = iv.uid;
			}
		}
	}

	for (auto iv : fake.entities) {
		if (iv.uid == A_snaffle) {
			double d = sqrt(euclid_dist2(vec2f(iv.x - A_dat->x, iv.y - A_dat->y)));
			if (d < 400) d = 400;
			fitness_value += c * (40 * 1000) / d;
		}
		if (iv.uid == B_snaffle) {
			double d = sqrt(euclid_dist2(vec2f(iv.x - B_dat->x, iv.y - B_dat->y)));
			if (d < 400) d = 400;
			fitness_value += c * (40 * 1000) / d;
		}
	}

	/*for (auto iv : fake.entities) {
		if(iv.type == entity::snaffle){
			float d_goal = sqrt(euclid_dist2(vec2f(iv.x - my_goal.first, iv.y - my_goal.second)));
			float d_A = sqrt(euclid_dist2(vec2f(iv.x - A_dat->x, iv.y - A_dat->y)));
			float d_B = sqrt(euclid_dist2(vec2f(iv.x - B_dat->x, iv.y - B_dat->y)));
			if (d_A < 400) d_A = 400;
			if (d_B < 400) d_B = 400;

			//if(d_A < d_B) fitness_value += c * (40 * d_goal) / (d_A);
			//else fitness_value += c * (40 * d_goal) / (d_B);
			fitness_value += c * (40 * d_goal) / (d_B + d_A);
			fitness_value += b * (40 * 1000) / d_goal;
		}
	}*/
	if (my_id == 0) fitness_value /= (fake.score_B + 1);
	else fitness_value /= (fake.score_A + 1);

	return fitness_value;
}

#define LOCAL 0

int main(){
	srand(time(0));

	if (!LOCAL) {
		int myTeamId; // if 0 you need to score on the right of the map, if 1 you need to score on the left
		std::cin >> myTeamId; std::cin.ignore();
		int my_mana = 0;
		int bludger_A_state(-1), bludger_B_state(-1);
		// game loop
		while (1) {
			
			gameState *state = new gameState;
			bool b__ = false;

			int entities; // number of entities still in game
			std::cin >> entities; std::cin.ignore();
			for (int i = 0; i < entities; i++) {
				std::string entityType;
				int pos_x, pos_y, vec_x, vec_y, entity_state(0), unique_id, parent_id(0);
				std::cin >> unique_id >> entityType >> pos_x >> pos_y >> vec_x >> vec_y >> entity_state; std::cin.ignore();
				
				entity::entity_type entity_type;
				if (myTeamId == 0) {
					if (entityType == "OPPONENT_WIZARD") parent_id = 1;
					else parent_id = 0;
				}
				else {
					if (entityType == "OPPONENT_WIZARD") parent_id = 0;
					else parent_id = 1;
				}
				

				if (entityType == "OPPONENT_WIZARD" || entityType == "WIZARD") entity_type = entity::wizard;
				else if (entityType == "SNAFFLE") entity_type = entity::snaffle;
				else {
					entity_type = entity::bludger;
					if(!b__){
						entity_state = bludger_A_state;
						b__ = true;
					}
					else {
						entity_state = bludger_B_state;
					}
				}
				if (entity_state == 1 && entity_type == entity::wizard) entity_state = 3;

				entity new_ent(entity{ unique_id, parent_id, entity_state, entity_type,
					(double)pos_x, (double)pos_y, (double)vec_x, (double)vec_y });

				std::cerr << unique_id << " " << entity_type << " " << pos_x << " " << pos_y << " " << vec_x << " " << vec_y
					<< " " << entity_state << " " << parent_id << std::endl;

				state->entities.push_back(new_ent);
			}

			entity new_ent{ -1, -1, -1, entity::immobile, 0, 1750, 0, 0 };
			state->entities.push_back(new_ent);
			new_ent = entity{ -1, -1, -1, entity::immobile, 0, 5750, 0, 0 };
			state->entities.push_back(new_ent);
			new_ent = entity{ -1, -1, -1, entity::immobile, 16000, 5750, 0, 0 };
			state->entities.push_back(new_ent);
			new_ent = entity{ -1, -1, -1, entity::immobile, 16000, 1750, 0, 0 };
			state->entities.push_back(new_ent);
			state->init();

			state->mana_A = my_mana;
			state->mana_B = my_mana;

			for (auto &iv : state->entities) {
				if (iv.type == entity::wizard && iv.state == 3) {
					for (auto kv : state->entities) {
						if (kv.type == entity::snaffle
							&& (kv.x - iv.x) < 0.1
							&& (kv.y - iv.y) < 0.1
							&& (kv.vx - iv.vx) < 0.1
							&& (kv.vy - iv.vy) < 0.1) {
							iv.grabbed_snaffle_uid = kv.uid;
						}
					}
				}
			}

			auto my_AI = new Genetic_A;
			my_AI->game_start(state, myTeamId);

			auto moves = my_AI->take_turn();

			if (moves.first.type == move::MOVE) {
				std::cout << "MOVE " << (int)std::round(moves.first.x) << " " <<
					(int)std::round(moves.first.y) << " " << 150 << std::endl;
			}
			else if (moves.first.type == move::THROW) {
				std::cout << "THROW " << (int)std::round(moves.first.x) << " " <<
					(int)std::round(moves.first.y) << " " << 500 << std::endl;
			}
			else if (moves.first.type == move::ACCIO) {
				std::cout << "ACCIO " << moves.first.t_uid << std::endl;
			}
			else if (moves.first.type == move::FLIPENDO) {
				std::cout << "FLIPENDO " << moves.first.t_uid << std::endl;
			}

			if (moves.second.type == move::MOVE) {
				std::cout << "MOVE " << (int)std::round(moves.second.x) << " " <<
					(int)std::round(moves.second.y) << " " << 150 << std::endl;
			}
			else if (moves.second.type == move::THROW) {
				std::cout << "THROW " << (int)std::round(moves.second.x) << " " <<
					(int)std::round(moves.second.y) << " " << 500 << std::endl;
			}
			else if (moves.second.type == move::ACCIO) {
				std::cout << "ACCIO " << moves.second.t_uid << std::endl;
			}
			else if (moves.second.type == move::FLIPENDO) {
				std::cout << "FLIPENDO " << moves.second.t_uid << std::endl;
			}

			dummy dumb;
			if (myTeamId == 0) state->play(moves, dumb.take_turn());
			else state->play(dumb.take_turn(), moves);

			if (myTeamId == 0) my_mana = state->mana_A;
			else my_mana = state->mana_B;

			b__ = false;
			for (auto iv : state->entities) {
				if (iv.type == entity::bludger) {
					if (!b__) {
						bludger_A_state = iv.state;
						b__ = true;
					}
					else {
						bludger_B_state = iv.state;
						break;
					}
				}
			}

			delete my_AI;
			delete state;
		}
	}
	else{
		int myTeamId = 0;
		gameState * state = new gameState;
		std::ifstream fin;
		fin.open("gamestates");
		int entity_num;

		fin >> entity_num;
		int unique_id, entity_type, pos_x, pos_y, vec_x, vec_y, entity_state, parent_id;
		for (int i(0); i < entity_num; ++i) {
			fin >> unique_id >> entity_type >> pos_x >> pos_y >> vec_x >> vec_y >> entity_state >> parent_id; fin.ignore();
			entity new_ent(entity{ unique_id, parent_id, entity_state, (entity::entity_type)entity_type,
				(double)pos_x, (double)pos_y, (double)vec_x, (double)vec_y });
			state->entities.push_back(new_ent);
		}

		for (auto &iv : state->entities) {
			if (iv.type == entity::wizard && iv.state == 3) {
				for (auto kv : state->entities) {
					if (kv.type == entity::snaffle
						&& (kv.x-iv.x) < 0.1
						&& (kv.y - iv.y) < 0.1
						&& (kv.vx - iv.vx) < 0.1
						&& (kv.vy - iv.vy) < 0.1) {
						iv.grabbed_snaffle_uid = kv.uid;
					}
				}
			}
		}

		entity new_ent{ -1, -1, -1, entity::immobile, 0, 1750, 0, 0 };
		state->entities.push_back(new_ent);
		new_ent = entity{ -1, -1, -1, entity::immobile, 0, 5750, 0, 0 };
		state->entities.push_back(new_ent);
		new_ent = entity{ -1, -1, -1, entity::immobile, 16000, 5750, 0, 0 };
		state->entities.push_back(new_ent);
		new_ent = entity{ -1, -1, -1, entity::immobile, 16000, 1750, 0, 0 };
		state->entities.push_back(new_ent);
		state->init();

		auto my_AI = new Genetic_A;
		my_AI->game_start(state, myTeamId);

		auto moves = my_AI->take_turn();

		if (moves.first.type == move::MOVE) {
			std::cout << "MOVE " << (int)std::round(moves.first.x) << " " <<
				(int)std::round(moves.first.y) << " " << 150 << std::endl;
		}
		else if (moves.first.type == move::THROW) {
			std::cout << "THROW " << (int)std::round(moves.first.x) << " " <<
				(int)std::round(moves.first.y) << " " << 500 << std::endl;
		}
		else if (moves.first.type == move::ACCIO) {
			std::cout << "ACCIO " << moves.first.t_uid << std::endl;
		}
		else if (moves.first.type == move::FLIPENDO) {
			std::cout << "FLIPENDO " << moves.first.t_uid << std::endl;
		}

		if (moves.second.type == move::MOVE) {
			std::cout << "MOVE " << (int)std::round(moves.second.x) << " " <<
				(int)std::round(moves.second.y) << " " << 150 << std::endl;
		}
		else if (moves.second.type == move::THROW) {
			std::cout << "THROW " << (int)std::round(moves.second.x) << " " <<
				(int)std::round(moves.second.y) << " " << 500 << std::endl;
		}
		else if (moves.second.type == move::ACCIO) {
			std::cout << "ACCIO " << moves.second.t_uid << std::endl;
		}
		else if (moves.second.type == move::FLIPENDO) {
			std::cout << "FLIPENDO " << moves.second.t_uid << std::endl;
		}

		std::cin >> entity_num;
	}
	return 0;
}