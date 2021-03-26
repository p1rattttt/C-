#include<list>
#include<iostream>
#include<memory>

template<size_t chunkSize = 1>
class FixedAllocator {
private:
    char* buffer_ = new char[1 << 27];
    size_t length = 0;
public:
    void* allocate(size_t need) {
        auto ptr = reinterpret_cast<void*>(buffer_ + length);
        length += need;
        return ptr;
    }

    void dealloc(void*) {
//        return nullptr;
    }

    ~FixedAllocator() {
        delete[] buffer_;
    }
};

template<typename T>
class FastAllocator {
private:
    static const size_t maxSize = 128;
    std::shared_ptr<FixedAllocator<>> allocator_;
public:
    using value_type = T;

    T* allocate(size_t need) {
        if (need * sizeof(T) <= maxSize) {
            return reinterpret_cast<T*>(allocator_->allocate(need * sizeof(T)));
        }
        return reinterpret_cast<T*>(::operator new(need * sizeof(T)));
    }

    void deallocate(T* ptr, size_t need) {
        if (need * sizeof(T) <= maxSize) {
            allocator_->dealloc(ptr);
            return;
        }
        ::operator delete(ptr);
    }

    std::shared_ptr<FixedAllocator<>> get_allocator() const {
        return allocator_;
    }

    FastAllocator(const FastAllocator& Other) {
        allocator_ = Other.allocator_;
    }

    template<typename U>
    FastAllocator(const FastAllocator<U>& Other) {
        allocator_ = Other.get_allocator();
    }

    FastAllocator() {
        allocator_ = std::make_shared<FixedAllocator<>>();
    }
};

template<typename T, typename U>
bool operator==(const FastAllocator<T>&, const FastAllocator<U>&) {
    return true;
}

template<typename T, typename U>
bool operator!=(const FastAllocator<T>&, const FastAllocator<U>&) {
    return false;
}

//template<typename T, typename Allocator = std::allocator<T>>
//using List = std::list<T, Allocator>;

template <typename T, typename Allocator = std::allocator<T>>
class List {
public:
    struct Node {
        Node* next;
        Node* prev;
        T x;
        explicit Node() : next(nullptr), prev(nullptr), x(T()) {}
        explicit Node(const T& x) : next(nullptr), prev(nullptr), x(x) {}
    };
    size_t length;
    Node* term;
    Allocator alloc;
    using AllocTraits = std::allocator_traits<Allocator>;
    using NodeAllocType = typename AllocTraits::template rebind_alloc<Node>;
    NodeAllocType NodeAlloc;
    using NodeAllocTraits = std::allocator_traits<NodeAllocType>;

    template<bool is_const>
    class TemplateIterator;
public:
    using iterator = TemplateIterator<false>;
    using const_iterator = TemplateIterator<true>;
    using reverse_iterator = std::reverse_iterator<iterator>;
    using const_reverse_iterator = std::reverse_iterator<const_iterator>;

    explicit List(const Allocator& alloc = Allocator()) : alloc(alloc) {
        length = 0;
        NodeAllocType newNodeAlloc;
        NodeAlloc = newNodeAlloc;
        term = NodeAllocTraits::allocate(NodeAlloc, 1);
        term->next = term;
        term->prev = term;
    }
    List(size_t count, const T& x, const Allocator& alloc = Allocator()) : alloc(alloc) {
        length = count;
        NodeAllocType newNodeAlloc;
        NodeAlloc = newNodeAlloc;
        term = NodeAllocTraits::allocate(NodeAlloc, 1);
        term->next = term;
        term->prev = term;
        for (size_t i = 0; i < count; ++i) {
            push_back(x);
        }
    }
    List(size_t count, const Allocator& alloc = Allocator()) : alloc(alloc) {
        length = count;
        NodeAllocType newNodeAlloc;
        NodeAlloc = newNodeAlloc;
        term = NodeAllocTraits::allocate(NodeAlloc, 1);
        term->next = term;
        term->prev = term;
        Node* cur = term;
        for (size_t i = 0; i < count; ++i) {
            Node* nw = NodeAllocTraits::allocate(NodeAlloc, 1);
            NodeAllocTraits::construct(NodeAlloc, nw);
            nw->prev = cur;
            cur->next = nw;
            nw->next = term;
            term->prev = nw;
            cur = nw;
        }
    }
    Allocator get_allocator() const {
        return alloc;
    }
    List(const List& x) {
        alloc = AllocTraits::select_on_container_copy_construction(x.alloc);
        NodeAlloc = NodeAllocTraits::select_on_container_copy_construction(x.NodeAlloc);
        Node* curOther = x.term->next;
        term = NodeAllocTraits::allocate(NodeAlloc, 1);
        term->next = term;
        term->prev = term;
        length = 0;
        for (size_t i = 0; i < x.length; ++i) {
            push_back(curOther->x);
            curOther = curOther->next;
        }
    }
    ~List() {
        Node* cur = term->next;
        while (cur != term) {
            Node* next = cur->next;
            NodeAllocTraits::destroy(NodeAlloc, cur);
            NodeAllocTraits::deallocate(NodeAlloc, cur, 1);
            cur = next;
        }
        NodeAllocTraits::deallocate(NodeAlloc, cur, 1);
    }
    List& operator=(const List& x) {
        if (this == &x)
            return *this;
        if (AllocTraits::propagate_on_container_copy_assignment::value && alloc != x.alloc) {
            alloc = x.alloc;
            NodeAlloc = x.NodeAlloc;
        }
        List nw = x;
        std::swap(nw.length, this->length);
        std::swap(nw.term, this->term);
        return *this;
    }
    void insert(const_iterator it, const T& x) {
        Node* iter = const_cast<Node*>(it.cur);
        Node* nw = NodeAllocTraits::allocate(NodeAlloc, 1);
        NodeAllocTraits::construct(NodeAlloc, nw, x);
        nw->prev = iter->prev;
        iter->prev->next = nw;
        nw->next = iter;
        iter->prev = nw;
        ++length;

    }
    void erase(const_iterator it) {
        Node* iter = const_cast<Node*>(it.cur);
        iter->next->prev = iter->prev;
        iter->prev->next = iter->next;
        --length;
        NodeAllocTraits::destroy(NodeAlloc, iter);
        NodeAllocTraits::deallocate(NodeAlloc, iter, 1);
    }

    size_t size() const {
        return length;
    }
    void push_back(const T& x) {
        Node* nw = NodeAllocTraits::allocate(NodeAlloc, 1);
        NodeAllocTraits::construct(NodeAlloc, nw, x);
        term->prev->next = nw;
        nw->next = term;
        nw->prev = term->prev;
        term->prev = nw;
        ++length;
    }
    void push_front(const T& x) {
        Node* nw = NodeAllocTraits::allocate(NodeAlloc, 1);
        NodeAllocTraits::construct(NodeAlloc, nw, x);
        term->next->prev = nw;
        nw->next = term->next;
        nw->prev = term;
        term->next = nw;
        ++length;
    }
    void pop_back() {
        Node* nw = term->prev;
        nw->prev->next = term;
        term->prev = nw->prev;
        NodeAllocTraits::destroy(NodeAlloc, nw);
        NodeAllocTraits::deallocate(NodeAlloc, nw, 1);
        --length;
    }
    void pop_front() {
        Node* nw = term->next;
        term->next = nw->next;
        nw->next->prev = term;
        NodeAllocTraits::destroy(NodeAlloc, nw);
        NodeAllocTraits::deallocate(NodeAlloc, nw, 1);
        --length;
    }
    iterator begin() const {
        return iterator(term->next);
    }
    iterator end() const {
        return iterator(term);
    }
    const_iterator cbegin() const {
        return const_iterator(term->next);
    }
    const_iterator cend() const {
        return const_iterator(term);
    }
    reverse_iterator rbegin() const {
        return reverse_iterator(term);
    }
    reverse_iterator rend() const {
        return reverse_iterator(term->next);
    }
    const_reverse_iterator crbegin() const {
        return const_reverse_iterator(term);
    }
    const_reverse_iterator crend() const {
        return const_reverse_iterator(term->next);
    }
};

template<typename T, typename Allocator>
template<bool is_const>
class List<T, Allocator>::TemplateIterator {
public:
    using type = List<T, Allocator>::Node*;
    type cur;
public:
    using iterator_category = std::bidirectional_iterator_tag;
    using value_type = T;
    using reference = std::conditional_t<is_const, const T&, T&>;
    using pointer = std::conditional_t<is_const, const T*, T*>;
    using difference_type = std::ptrdiff_t;

    TemplateIterator(Node* cur) : cur(cur) {}

    operator TemplateIterator<1>() {
        return TemplateIterator<1>(cur);
    }
    bool operator==(const TemplateIterator& x) const {
        return cur == x.cur;
    }
    bool operator!=(const TemplateIterator& x) const {
        return cur != x.cur;
    }
    std::conditional_t<is_const, const T&, T&> operator*() const {
        return cur->x;
    }
    std::conditional_t<is_const, const T*, T*> operator->() const {
        return &(cur->x);
    }
    TemplateIterator& operator++() {
        cur = cur->next;
        return *this;
    }
    TemplateIterator& operator--() {
        cur = cur->prev;
        return *this;
    }
    TemplateIterator operator++(int) {
        TemplateIterator copy = *this;
        cur = cur->next;
        return copy;
    }
    TemplateIterator operator--(int) {
        TemplateIterator copy = *this;
        cur = cur->prev;
        return copy;
    }
};
