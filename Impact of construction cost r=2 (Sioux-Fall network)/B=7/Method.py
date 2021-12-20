from Data import Read_data
from gurobipy import *
import copy
class Solve:
    def __init__(self):
        #input
        self.multiplier=0
        data=Read_data(self.multiplier)
        self.node_list, self.link_list, self.candidate_link_list, self.OD_pair_list, \
        self.g_number_of_links, self.g_number_of_nodes, self.g_number_of_candidate_links, self.g_number_of_ODs=data.read_candidate_links()
        self.reliability=2
        self.construction_budget=7
        self.construction_cost=1
        self.iteration_times=20
        self.acceptable_gap=0.000

        #output
        # self.record_multiplier_miu=[] #optimal dual variable
        self.local_LB = []
        self.local_UB = []
        self.global_LB = []
        self.global_UB = []
        self.solutions_of_routing_subproblem=[]
        self.solutions_of_RMP=[]
        self.obtained_optimal_dual_value=[]

    def g_solving_RNDP_by_BD(self):
        print("Solving...")
        self.g_define_the_RMP()
        for i in range(self.iteration_times):
            # self.record_multiplier_miu.append([])
            self.solutions_of_routing_subproblem.append([])
            self.obtained_optimal_dual_value.append([])
            self.local_LB.append(0)
            self.local_UB .append(0)
            self.global_LB .append(-100000)
            self.global_UB .append(100000)
            print("Iteration: {}".format(i+1))

#I. solve the RMP
            print("RMP_{}".format(i+1))
            self.construction=[]
            self.RMP.optimize()
            values=self.RMP.getVars()
            for candidate_link_index in range(self.g_number_of_candidate_links):
                result=round(values[candidate_link_index].x)
                self.construction.append(result)
            self.global_LB[i]=self.RMP.objval
            self.local_LB[i]=self.RMP.objval
            # print(self.RMP.objval)
            self.solutions_of_RMP.append(self.construction)

#II. solve the routing subproblems according to the construction
            for od_index in range(self.g_number_of_ODs):
                od_pair = self.OD_pair_list[od_index]
                #solve this problem by LR and output a BD cut
                node_seq, candidate_travel_flag, global_lb, global_ub,optimal_dual_variable=self.g_solving_RSP(od_index,od_pair,i)
                self.local_UB[i]+=global_ub
                self.solutions_of_routing_subproblem[i].append(node_seq)
                self.obtained_optimal_dual_value[i].append(optimal_dual_variable)
            # print(self.local_UB[i])

            if i==0:
                self.global_UB[0]=self.local_UB[0]
            else:
                self.global_UB[i] = min(self.local_UB[i], self.global_UB[i-1])

            self.RMP.update()

#III. Terminal condition test
            gap=(self.global_UB[i]-self.global_LB[i])/self.global_UB[i]
            if gap<=self.acceptable_gap or i==self.iteration_times-1:
                self.consumed_iteration=i+1
                self.RMP.write("RMP.lp")
                break



    def g_define_the_RMP(self):
        self.RMP=Model()
        self.RMP.setParam('OutputFlag', 0)
        #1.add variable
        for link in self.candidate_link_list:
            from_node_id=link.from_node_id
            to_node_id=link.to_node_id
            self.RMP.addVar(vtype=GRB.BINARY,name="y_{}_{}".format(from_node_id,to_node_id))
        for od_index in range(self.g_number_of_ODs):
            self.RMP.addVar(lb=0,obj=1,vtype=GRB.CONTINUOUS,name="z_{}".format(od_index))
        self.RMP.update()

        #2.add constraints
        #construction budget
        expr=LinExpr()
        for link in self.candidate_link_list:
            from_node_id=link.from_node_id
            to_node_id=link.to_node_id
            name=self.RMP.getVarByName("y_{}_{}".format(from_node_id,to_node_id))
            expr.addTerms(1,name)
        # self.RMP.addConstr(expr,GRB.LESS_EQUAL,self.construction_budget,name="budget")
        self.RMP.addConstr(expr, GRB.LESS_EQUAL, self.construction_budget, name="budget")

        self.RMP.write("RMP.lp")
        print("RMP is initialized!!!")



    def g_solving_RSP(self,od_index,od_pair,i):
        #input:network,od,pair
        #output:global lower bound,the most reliabile shortest path [node sequence] and a Benders cut
        #1.Initilization;find the interval of the variance
        global_lb=[]
        global_ub=[]
        number_of_iterations=30
        accsptable_gap=0.005
        multiplier_gama=0
        optimal_solution=None
        multiplier_list=[]
        #record the information of BD cut when the optimal lb is obtained
        optimal_multiplier=None
        optimal_L2=None
        optimal_dual_variable_of_L1=None

        values,obj,dual_variable_list=self.g_solving_SP(od_index,od_pair,multiplier_gama,"1")

        max_variance = 0
        for link in self.link_list:
            variance = link.travel_time_variance
            index = link.link_id
            if round(values[index].x) == 1:
                max_variance += variance

        #2.solve two subproblems
        for j in range(number_of_iterations):
            multiplier_list.append(multiplier_gama)
            # print("RSP_{}".format(j+1))
            #Subproblem I
            values,obj,dual_variable_list=self.g_solving_SP(od_index,od_pair,multiplier_gama,"2")

            path_mean=0
            path_variance=0
            path_multiplier=0
            for link in self.link_list:
                mean=link.travel_time_mean
                variance = link.travel_time_variance
                index = link.link_id
                if round(values[index].x) == 1:
                    path_variance += variance
                    path_mean+=mean
                    if link.link_type==1:#candidate link
                        path_multiplier+=link.base_profit_for_lagrangian[od_index]
            local_lb=obj

            # Subproblem II
            L1=0
            L2=self.reliability*(max_variance)**0.5-multiplier_gama*max_variance
            if L1<L2:
                local_lb+=L1
                z_value=0
            else:
                local_lb+=L2
                z_value=max_variance

            #3.update bounds
            #generate a feasible solution for RSP
            local_ub = path_mean + self.reliability * path_variance ** 0.5

            if j==0:
                global_lb.append(local_lb)
                global_ub.append(local_ub)
                optimal_solution=values
                optimal_multiplier=multiplier_gama
                optimal_obj=obj
                optimal_L2=min(L1,L2)
                optimal_dual_variable_of_L1=dual_variable_list
            else:
                if local_lb>global_lb[-1]:
                    global_lb.append(local_lb)
                    optimal_multiplier = multiplier_gama
                    optimal_obj = obj
                    optimal_L2 = min(L1, L2)
                    optimal_dual_variable_of_L1 = dual_variable_list
                else:
                    global_lb.append(global_lb[-1])

                if local_ub<global_ub[-1]:
                    global_ub.append(local_ub)
                    optimal_solution=values
                else:
                    global_ub.append(global_ub[-1])

            #4.update multipliers
            if path_variance!=z_value:
                multiplier_gama+=(global_ub[-1]-local_lb)/(path_variance-z_value)

            #5.terminal conditions and add a Benders optimality cut
            terminal_flag=0
            if global_ub[-1]!=0:
                gap=(global_ub[-1]-global_lb[-1])/global_ub[-1]
                if gap<accsptable_gap:
                    terminal_flag=1

            if j == number_of_iterations - 1:
                terminal_flag = 1
                # print("warning")

            if terminal_flag==1:
                node_seq, candidate_travel_flag = self.values_transition(values, od_pair)
                # add a Benders cut
                #z_k
                expr=LinExpr()
                name=self.RMP.getVarByName(name="z_{}".format(od_index))
                expr.addTerms(1,name)
                #f_i_j
                for link in self.candidate_link_list:
                    from_node_id = link.from_node_id
                    to_node_id = link.to_node_id
                    link_index=self.candidate_link_list.index(link)
                    name=self.RMP.getVarByName(name="y_{}_{}".format(from_node_id, to_node_id))
                    value=optimal_dual_variable_of_L1[link_index]
                    expr.addTerms(value,name)
                rhs=optimal_L2+optimal_obj#optimal_dual_variable_of_L1[-1]-optimal_dual_variable_of_L1[-2]
                self.RMP.addLConstr(expr,GRB.GREATER_EQUAL,rhs,name="CUT_{}_{}".format(i,od_index))

                optimal_dual_variable_of_L1.append(optimal_multiplier)
                # print(j+1)
                # print(multiplier_gama)
                # print(global_ub[-1], global_lb[-1])
                return node_seq, candidate_travel_flag, global_lb[-1], global_ub[-1],optimal_dual_variable_of_L1


    def g_solving_SP(self,od_index,od_pair,multiplier_gama,Flag):
        # least expected path
        self.sp = Model("sp")
        self.sp.setParam('OutputFlag', 0)
        expr = LinExpr()
        for link in self.link_list:
            name = "x_{}_{}".format(link.from_node_id, link.to_node_id)
            name = self.sp.addVar(vtype=GRB.CONTINUOUS, name=name, lb=0, ub=1)
            value=0
            if Flag == "1": #and link.construction_Flag==1:
                value += link.travel_time_mean
            if Flag == "2":#and link.construction_Flag==1:
                value += link.travel_time_mean+link.travel_time_variance * multiplier_gama
            #when not be constrcuted, a large cost is added
            # if Flag == "2_1" and link.construction_Flag==0:
            #     value += 10000
            # if Flag == "2_2" and link.construction_Flag==0:
            #     value += 10000

            expr.addTerms(value, name)

        self.sp.setObjective(expr, GRB.MINIMIZE)
        self.sp.update()

        #add construction constrainsts
        for link_index in range(self.g_number_of_candidate_links):
            link=self.candidate_link_list[link_index]
            name=self.sp.getVarByName("x_{}_{}".format(link.from_node_id, link.to_node_id))
            value = self.construction[link_index]
            self.sp.addConstr(name, GRB.LESS_EQUAL, value,name="const_{}_{}".format(link.from_node_id, link.to_node_id))

        # Flow balance
        for node in self.node_list:
            expr = LinExpr()

            # Flow out
            for outbound_link in node.outbound_links_list:
                name = self.sp.getVarByName("x_{}_{}".format(outbound_link.from_node_id, outbound_link.to_node_id))
                expr.addTerms(1, name)

            # Flow in
            for inbound_link in node.inbound_links_list:
                name = self.sp.getVarByName("x_{}_{}".format(inbound_link.from_node_id, inbound_link.to_node_id))
                expr.addTerms(-1, name)

            if node.node_id ==od_pair[0]:
                self.sp.addConstr(expr, GRB.EQUAL, 1, name="Node_{}".format(od_pair[0]))

            if node.node_id ==od_pair[1]:
                self.sp.addConstr(expr, GRB.EQUAL, -1, name="Node_{}".format(od_pair[1]))

            if node.node_id != od_pair[0] and node.node_id != od_pair[1]:
                self.sp.addConstr(expr, GRB.EQUAL, 0, name="Node_{}".format(node.node_id))

        # self.sp.write("sp.lp")
        self.sp.optimize()
        values=self.sp.getVars()
        obj=self.sp.objval
        self.sp.update()

        #we should also return the best dual values
        dual_variable_list=[] #od,candidate links
        #candidate links
        for link_index in range(self.g_number_of_candidate_links):
            link=self.candidate_link_list[link_index]
            const = self.sp.getConstrByName(name="const_{}_{}".format(link.from_node_id, link.to_node_id))
            pi=round(const.pi,3)*-1
            dual_variable_list.append(pi)
        # origin
        const = self.sp.getConstrByName(name="Node_{}".format(od_pair[0]))
        pi = round(const.pi, 3) * -1
        dual_variable_list.append(pi)
        # destination
        const = self.sp.getConstrByName(name="Node_{}".format(od_pair[1]))
        pi = round(const.pi, 3) * -1
        dual_variable_list.append(pi)

        return values, obj,dual_variable_list

    def values_transition(self,values,od_pair):
        #input: values;output:node seq,use seq
        # #path
        path_links = {}
        candidate_travel_flag=[]
        for link in self.link_list:
            link_index = link.link_id
            from_node = link.from_node_id
            to_node = link.to_node_id
            if round(values[link_index].x) == 1:
                path_links[from_node] = to_node
            link_type = link.link_type

            if link_type==1 and round(values[link_index].x) == 1:
                candidate_travel_flag.append(1)
            if link_type==1 and round(values[link_index].x) == 0:
                candidate_travel_flag.append(0)

        node_seq = [od_pair[0]+1]
        current_node = od_pair[0]
        while current_node != od_pair[1]:
            current_node = path_links[current_node]
            node_seq.append(current_node+1)
        return node_seq,candidate_travel_flag



    def output_results(self,spend_time):
        with open("output_gap.csv","w") as fl:
            fl.write("iteration,local_LB,local_UB,LB,UB,gap\n")
            for i in range(len(self.global_UB)):
                local_LB=round(self.local_LB[i],3)
                local_UB=round(self.local_UB[i],3)
                LB=round(self.global_LB[i],3)
                UB=round(self.global_UB[i],3)
                gap=0
                if UB!=0:
                    gap=round((UB-LB)/UB,3)
                fl.write(str(i+1)+","+str(local_LB)+","+str(local_UB)+","+str(LB)+","+str(UB)+","+str(gap)+"\n")
            fl.write("CPU time: {}".format(spend_time))

        with open("output_solution_of_routing_subproblems.txt","w") as fl:
            fl.write("iteration,OD,path\n")
            for i in range(len(self.solutions_of_routing_subproblem)):
                solutions=self.solutions_of_routing_subproblem[i]
                for k in range(self.g_number_of_ODs):
                    od_pair=self.OD_pair_list[k]
                    fl.write(str(i+1)+","+"{}_{},".format(od_pair[0]+1,od_pair[1]+1))
                    fl.write(str(solutions[k]))
                    fl.write("\n")

        with open("outout_solution_of_KS.txt","w") as fl:
            fl.write("iteration,solution\n")
            for i in range(len(self.global_UB)):
                solution=[]
                result=self.solutions_of_RMP[i]
                for index in range(self.g_number_of_candidate_links):
                    if result[index]==1:
                        link=self.candidate_link_list[index]
                        from_node=link.from_node_id+1
                        to_node=link.to_node_id+1
                        solution.append((from_node,to_node))

                fl.write(str(i+1)+",")
                fl.write(str(solution)+"\n")

        with open("output_the_optimal_dual.csv","w") as fl:
            fl.write("iteration,OD,")
            for candidate_link in range(self.g_number_of_candidate_links):
                fl.write(str(candidate_link)+",")
            fl.write("destination,origin,multiplier_gama")
            fl.write("\n")

            for i in range(len(self.obtained_optimal_dual_value)):
                solutions = self.obtained_optimal_dual_value[i]
                for k in range(self.g_number_of_ODs):
                    od_pair = self.OD_pair_list[k]
                    fl.write(str(i + 1) + "," + "{}_{},".format(od_pair[0]+1, od_pair[1]+1))

                    for value in solutions[k]:
                        value=round(value,3)
                        fl.write(str(value)+",")
                    fl.write("\n")


        with open("output_the_optimal_dual_version_II.csv","w") as fl:
            fl.write("iteration,OD,")
            for candidate_link in range(self.g_number_of_candidate_links):
                fl.write(str(candidate_link)+",")
            fl.write("destination,origin,multiplier_gama")
            fl.write("\n")

            for k in range(self.g_number_of_ODs):
                for i in range(len(self.obtained_optimal_dual_value)):
                    solutions = self.obtained_optimal_dual_value[i]
                    od_pair = self.OD_pair_list[k]
                    fl.write(str(i + 1) + "," + "{}_{},".format(od_pair[0]+1, od_pair[1]+1))

                    for value in solutions[k]:
                        value=round(value,3)
                        fl.write(str(value)+",")

                    fl.write("\n")