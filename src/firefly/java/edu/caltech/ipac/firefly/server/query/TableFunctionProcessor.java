package edu.caltech.ipac.firefly.server.query;

import edu.caltech.ipac.firefly.data.FileInfo;
import edu.caltech.ipac.firefly.data.Param;
import edu.caltech.ipac.firefly.data.ServerParams;
import edu.caltech.ipac.firefly.data.TableServerRequest;
import edu.caltech.ipac.firefly.server.db.DbAdapter;
import edu.caltech.ipac.firefly.server.db.DbInstance;
import edu.caltech.ipac.firefly.server.db.EmbeddedDbUtil;
import edu.caltech.ipac.firefly.server.db.spring.JdbcFactory;
import edu.caltech.ipac.firefly.server.util.QueryUtil;
import edu.caltech.ipac.firefly.server.util.ipactable.DataGroupPart;
import edu.caltech.ipac.firefly.server.util.ipactable.TableDef;
import edu.caltech.ipac.util.DataGroup;
import edu.caltech.ipac.util.StringUtils;
import org.apache.commons.codec.digest.DigestUtils;

import java.io.File;
import java.util.SortedSet;
import java.util.TreeSet;

import static edu.caltech.ipac.firefly.data.TableServerRequest.FILTERS;

public abstract class TableFunctionProcessor extends EmbeddedDbProcessor {
    public static final String SEARCH_REQUEST = "searchRequest";


    abstract protected DataGroup fetchData(TableServerRequest treq, File dbFile, DbAdapter dbAdapter) throws DataAccessException;

    /**
     *  recreate the table database if does not exists.  otherwise, return the table's database
     */
    @Override
    public FileInfo createDbFile(TableServerRequest treq) throws DataAccessException {
        TableServerRequest sreq = getSearchRequest(treq);
        File dbFile = EmbeddedDbUtil.getDbFile(sreq);
        if (dbFile == null || !dbFile.canRead()) {
            sreq.setPageSize(1);
            sreq.setStartIndex(0);
            new SearchManager().getDataGroup(sreq).getData();
        }
        EmbeddedDbUtil.setDbMetaInfo(treq, DbAdapter.getAdapter(treq), dbFile);
        return new FileInfo(dbFile);
    }

    /**
     * generate stats for the given search request if not exists.  otherwise, return the stats
     */
    @Override
    protected DataGroupPart getResultSet(TableServerRequest treq, File dbFile) throws DataAccessException {

        String resTblName = getResultSetTable(treq);

        DbAdapter dbAdapter = DbAdapter.getAdapter(treq);
        DbInstance dbInstance =  dbAdapter.getDbInstance(dbFile);
        String tblExists = String.format("select count(*) from %s", resTblName);
        try {
            JdbcFactory.getSimpleTemplate(dbInstance).queryForInt(tblExists);
        } catch (Exception e) {
            // does not exists.. fetch data and populate
            DataGroup data = fetchData(treq, dbFile, dbAdapter);
            EmbeddedDbUtil.createDataTbl(dbFile, data, dbAdapter, resTblName);
            EmbeddedDbUtil.createDDTbl(dbFile, data, dbAdapter, resTblName);
            EmbeddedDbUtil.createMetaTbl(dbFile, data, dbAdapter, resTblName);
        }
        treq.setParam(TableServerRequest.SQL_FROM, resTblName);
        String sql = String.format("%s %s %s", dbAdapter.selectPart(treq), dbAdapter.fromPart(treq), dbAdapter.wherePart(treq));
        sql = dbAdapter.translateSql(sql);

        DataGroup dg = EmbeddedDbUtil.runQuery(dbAdapter, dbFile, sql, resTblName);
        TableDef tm = new TableDef();
        tm.setStatus(DataGroupPart.State.COMPLETED);
        return new DataGroupPart(tm, dg, treq.getStartIndex(), dg.size());
    }

    protected TableServerRequest getSearchRequest(TableServerRequest treq) throws DataAccessException {
        String searchRequestJson = treq.getParam(SEARCH_REQUEST);
        if (searchRequestJson == null) {
            throw new DataAccessException("Action failed: " + SEARCH_REQUEST + " is missing");
        }
        TableServerRequest sreq = QueryUtil.convertToServerRequest(searchRequestJson);
        if (sreq.getRequestId() == null) {
            throw new DataAccessException("Action failed: " + SEARCH_REQUEST + " must contain " + ServerParams.ID);
        }
        return sreq;
    }

    /**
     * returns the table name of the resultset.  the same request should return the same table name so that it
     * does not need to be recreated.
     * @param treq
     * @return
     * @throws DataAccessException
     */
    protected String getResultSetTable(TableServerRequest treq) throws DataAccessException {
        TableServerRequest sreq = getSearchRequest(treq);
        TreeSet<Param> params = new TreeSet<>();
        if (sreq.getFilters() != null && sreq.getFilters().size() > 0) {
            params.add(new Param(FILTERS, TableServerRequest.toFilterStr(sreq.getFilters())));
        }
        params.addAll(treq.getSearchParams());
        String id = StringUtils.toString(params, "|");
        return String.format("%s_data_%s", getResultSetTablePrefix(), DigestUtils.md5Hex(id));
    }

    abstract protected String getResultSetTablePrefix();

}










