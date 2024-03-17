#include "JobQueue.h"


bool JobAbstractQueue::compareNode(JobData *a, JobData *b) {
    return a->node_data->node->ub > b->node_data->node->ub;
}

JobData *JobAbstractQueue::pop() {
    JobData *node = queue.front();
    queue.pop_front();
    return node;
}

bool JobAbstractQueue::empty() {
    return queue.empty();
}

void JobAbstractQueue::sort() {
    std::sort(queue.begin(), queue.end(), compareNode);
}

Node *JobAbstractQueue::getMinUb() {
    double min_ub = std::numeric_limits<double>::infinity();
    Node *min_node = nullptr;
    for (auto current : queue) {
        if (current->node_data->node->ub < min_ub) {
            min_ub = current->node_data->node->ub;
            min_node = current->node_data->node;
        }
    }
    return min_node;
}

Node *JobAbstractQueue::getMaxUb() {
    double max_ub = - std::numeric_limits<double>::infinity();
    Node *max_node = nullptr;
    for (auto current : queue) {
        if (current->node_data->node->ub > max_ub) {
            max_ub = current->node_data->node->ub;
            max_node = current->node_data->node;
        }
    }
    return max_node;
}

void JobAbstractQueue::print() {
    for (auto &elem : queue) {
        std::cout << elem->node_data->node->ub << " ";
    }
    std::cout << "\n";
}

int JobAbstractQueue::getSize() {
    return (int) queue.size();
}

void JobQueue::push(JobData *node) {
    queue.push_back(node);
}

void JobStack::push(JobData *node) {
    queue.push_front(node);
}

void JobPriorityQueue::push(JobData *node) {
    queue.push_back(node);
    sort();
}
