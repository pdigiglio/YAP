#include "DataSet.h"

#include "DataPoint.h"
#include "Exceptions.h"
#include "Model.h"

namespace yap {

//-------------------------
DataSet::DataSet(const Model& m) :
    DataPartitionBlock(m.dataAccessors()),
    Model_(&m)
{
}

//-------------------------
DataSet::DataSet(const DataSet& other) :
    DataPartitionBlock(other),
    DataPoints_(other.DataPoints_),
    Model_(other.Model_)
{
    assertDataPointOwnership();
}

//-------------------------
DataSet::DataSet(DataSet&& other) :
    DataPartitionBlock(std::move(other)),
    DataPoints_(std::move(other.DataPoints_)),
    Model_(std::move(other.Model_))
{
    assertDataPointOwnership();
}

//-------------------------
DataSet& DataSet::operator=(const DataSet& other)
{
    DataPartitionBlock::operator=(other);
    Model_ = other.Model_;
    DataPoints_ = other.DataPoints_;
    assertDataPointOwnership();
    return *this;
}

//-------------------------
DataSet& DataSet::operator=(DataSet&& other)
{
    DataPartitionBlock::operator=(std::move(other));
    Model_ = std::move(other.Model_);
    DataPoints_ = std::move(other.DataPoints_);
    assertDataPointOwnership();
    return *this;
}

//-------------------------
void DataSet::swap(DataSet& other)
{
    using std::swap;
    swap(static_cast<DataPartitionBlock&>(*this), static_cast<DataPartitionBlock&>(other));
    swap(Model_, other.Model_);
    swap(DataPoints_, other.DataPoints_);
    assertDataPointOwnership();
    other.assertDataPointOwnership();
}

//-------------------------
void DataSet::assertDataPointOwnership()
{
    for (auto& d : DataPoints_)
        d.DataSet_ = this;
}

//-------------------------
bool DataSet::consistent(const DataPoint& d) const
{
    return points().empty() or equalStructure(points().front(), d);
}

//-------------------------
void DataSet::addEmptyPoint()
{
    if (!model())
        throw exceptions::Exception("Model unset or deleted", "DataSet::add");

    DataPoints_.emplace_back(*this);
    auto& d = DataPoints_.back();

    if (!consistent(d))
        throw exceptions::Exception("produced inconsistent data point", "Model::addDataPoint");
}

//-------------------------
void DataSet::addEmptyPoints(size_t n)
{
    for (size_t i = 0; i < n; ++i)
        addEmptyPoint();
}

//-------------------------
void DataSet::add(const std::vector<FourVector<double> >& P)
{
    addEmptyPoint();
    DataPoints_.back().setFinalStateMomenta(P);
}

}
