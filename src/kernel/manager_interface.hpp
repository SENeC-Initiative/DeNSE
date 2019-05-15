#ifndef MANAGER_INTERFACE_H
#define MANAGER_INTERFACE_H

namespace growth
{

/**
 * Interface for kernel manager classes.
 *
 * This class defines the common interface for all manager classes
 *
 * @note Each manager shall be instantiated only once. Therefore, copy
 * constructor and assignment operator are declared private and not implemented.
 *
 * @ingroup KernelManagers
 */
class ManagerInterface
{
  private:
    ManagerInterface(ManagerInterface const &); // do not implement
    void operator=(ManagerInterface const &);   // do not implement

  public:
    ManagerInterface() {}

    virtual ~ManagerInterface() {}

    /*
     * Prepare manager for operation.
     *
     * After this method has completed, the manager should be completely
     * initialized and "ready for action".
     *
     * @note Initialization of any given manager may depend on other
     * managers having been initialized before. KernelManager::initialize()
     * is responsible for calling the initialization routines on the
     * specific managers in correct order.
     *
     * @see finalize()
     */
    virtual void initialize() = 0;

    /*
     * Take down manager after operation.
     *
     * After this method has completed, all dynamic data structures created
     * by the manager shall be deallocated and containers emptied. Plain
     * variables need not be reset.
     *
     * @note Finalization of any given manager may depend on other
     * managers not having been finalized yet. KernelManager::finalize()
     * is responsible for calling the initialization routines on the
     * specific managers in correct order, i.e., the opposite order of
     * initialize() calls.
     *
     * @see initialize()
     */
    virtual void finalize() = 0;
};
}

#endif /* MANAGER_INTERFACE_H */
