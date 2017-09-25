/*
 * License information at https://github.com/Caltech-IPAC/firefly/blob/master/License.txt
 */
package edu.caltech.ipac.firefly.server.db;

import edu.caltech.ipac.firefly.data.TableServerRequest;

import java.io.File;

/**
 * @author loi
 * @version $Id: DbInstance.java,v 1.3 2012/03/15 20:35:40 loi Exp $
 */
public class H2DbAdapter extends BaseDbAdapter{

    public String getName() {
        return H2;
    }

    public DbInstance getDbInstance(File dbFile) {
        String dbUrl = String.format("jdbc:h2:%s;CACHE_SIZE=1048576;LOG=0;UNDO_LOG=0;MVCC=true", dbFile.getPath());
        return new EmbeddedDbInstance(getName(), dbUrl, "org.h2.Driver");
    }

    public File getStorageFile(File dbFile) {
        return dbFile == null ? null : new File(dbFile.getParent(), dbFile.getName() + ".mv.db");
    }

    public String createTableFromSelect(String tblName, String selectSql) {
        return String.format("CREATE TABLE IF NOT EXISTS %s AS (%s)", tblName, selectSql);
    }

}
